#!/usr/bin/env nextflow
nextflow.preview.output = true
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

log.info """
    |            #################################################
    |            #    _  _                             _         #
    |            #   | \\| |  ___  __ __  ___   _ __   (_)  __    #
    |            #   | .` | / -_) \\ \\ / / _ \\ | '  \\  | | (_-<   #
    |            #   |_|\\_| \\___| /_\\_\\ \\___/ |_|_|_| |_| /__/   #
    |            #                                               #
    |            #################################################
    |
    | rna-preprocessing: Estimate the expression level of features (genes or transcripts) and perform the associated QC. Starting from cleaned reads and one pre-built reference database in standard format (cf. Readme.md).
    | 
    |""".stripMargin()

if (params.help) {
  log.info paramsHelp("nextflow run nexomis/rna-preprocessing --input </path/to/samplesheet> --kallisto_idx </path/to/kallisto/txome/file.index> [args]")
  exit 0
}

validateParameters()
log.info paramsSummaryLog(workflow)

file(params.out_dir + "/nextflow").mkdirs()

def parse_sample_entry(it) {
  def type = "SR"
  def r1_file = file(it[1], checkIfExists: true)
  def files = [r1_file]
  if (it[2] && !it[2].isEmpty() ) {
    def r2_file = file(it[2], checkIfExists: true)
    files << r2_file
    type = "PE"
  }
  if (r1_file.toString().toLowerCase().endsWith("sfq")) {
    type = "sfq"
  }
  def meta = ["id": it[0], "kallisto_idx": it[3],
      "strand": it[4], "is_3prime": it[5],
      "frag_size": it[6],"frag_size_sd": it[7], "read_type": type
    ]
  return [meta, files]
}

def parse_kallisto_idx(it) {
  def is_index = false
  def is_fasta = false
  def ref_file
  if (it[1] && !it[1].isEmpty()) {
    ref_file = file(it[1], checkIfExists: true)
    is_index = true
  } else if (it[2] && !it[2].isEmpty()) {
    ref_file = file(it[2], checkIfExists: true)
    is_fasta = true
  }

  def meta = ["id": it[0], "is_index": is_index, "is_fasta": is_fasta]
  return [meta, ref_file]
}

def make_kallisto_args(meta) {
  def kallisto_args = []

  if (meta.strand == "fr-stranded") {
    kallisto_args << "--fr-stranded"
  } else if (meta.strand == "rf-stranded") {
    kallisto_args << "--rf-stranded"
  }
  if (meta.read_type == "SR") {
    kallisto_args << "--single"
  }
  if (meta.is_3prime) {
    kallisto_args << "--single-overhang"
  }
  def has_frag_size = false
  if (meta.frag_size > 0) {
    kallisto_args << "--fragment-length " + meta.frag_size
    if (meta.frag_size_sd > 0) {
      kallisto_args << "--sd "+ meta.frag_size_sd
      has_frag_size = true
    }
  }

  if ((! has_frag_size) && (meta.read_type == "SR")) {
    if (meta.is_3prime) {
      // add fake value for 3 prime it is not set
      kallisto_args << "--fragment-length 250"
      kallisto_args << "--sd 50"
    } else {
      error "Fragment size and standard deviation are required for single-end reads in sample ${meta.id}. Please provide frag_size and frag_size_sd in the samplesheet."
    }
  }
  
  return kallisto_args.join(" ")
}

include { PRIMARY } from './modules/subworkflows/primary/main.nf'
include { RNA_PREPROCESSING } from './modules/subworkflows/rna_preprocessing/main.nf'
include { KALLISTO_INDEX } from './modules/process/kallisto/index/main.nf'

workflow {
  // START PARSING SAMPLE SHEET
  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | set { rawInputs }

  // START PRIMARY
  if (params.skip_primary) {
    trimmedInputs = rawInputs
  } else {
    if ( params.kraken2_db == null ) {
      error "kraken2_db argument required for primary analysis"
    }

    Channel.fromPath(params.kraken2_db, type: "dir", checkIfExists: true)
    | map {[["id": "kraken_db"], it]}
    | collect
    | set {dbPathKraken2}
    
    taxDir = Channel.fromPath(params.tax_dir, type: 'dir', checkIfExists: true)
    numReads = Channel.value(params.num_reads_sample_qc)

    PRIMARY(rawInputs, dbPathKraken2, taxDir, numReads)
    
    trimmedInputs = PRIMARY.out.trimmed
  }
  // END PRIMARY

  // START RNA PREPROCESSING
  trimmedInputs
  | map {
    it[0].kallisto_args = make_kallisto_args(it[0])
    return it
  }
  | set { inputsForKallisto }

  Channel.fromSamplesheet("kallisto_idx")
  | map {
      return parse_kallisto_idx(it)
    }
  | set { kallistoIndex }

  kallistoIndex
  | branch {
    build: it[0].is_fasta
    existing: it[0].is_index
  }
  | set { kallistoIndexBranched }

  KALLISTO_INDEX(kallistoIndexBranched.build)
  
  kallistoIndexFinal = KALLISTO_INDEX.out.idx.concat(kallistoIndexBranched.existing)

  multiqcYml = Channel.fromPath(projectDir + "/files/align_multiqc.yml")

  RNA_PREPROCESSING(inputsForKallisto, kallistoIndexFinal, multiqcYml)

  publish:
  PRIMARY.out.trimmed                         >> 'trimmed_and_filtered'
  PRIMARY.out.fastqc_trim_html                >> 'fastqc_trim'
  PRIMARY.out.fastqc_raw_html                 >> 'fastqc_raw'
  PRIMARY.out.multiqc_html                    >> 'multiqc'
  RNA_PREPROCESSING.out.kallisto_h5           >> 'kallisto'
  RNA_PREPROCESSING.out.multiqc               >> 'align_multiqc'
}

output {
  'trimmed_and_filtered' {
    enabled params.save_fastp
  }
  'fastqc_trim' {
    path 'qc/fastqc/trimmed'
  }
  'fastqc_raw' {
    path 'qc/fastqc/raw'
  }
  'multiqc' {
    path 'qc/multiqc/primary'
  }
  'kallisto' {
    path 'kallisto'
  }
  'align_multiqc' {
    path 'qc/multiqc/align'
  }
}
