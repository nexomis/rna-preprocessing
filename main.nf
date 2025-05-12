#!/usr/bin/env nextflow
nextflow.preview.output = true
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'



def parse_sample_entry(it) {
    // Define the order of fields as per assets/input_schema.json
    def field_names = [
        "sample_name", "path_r1", "path_r2", "reference", "strand",
        "is_3prime", "frag_size", "frag_size_sd", "r1_adapter",
        "r2_adapter", "trim_poly_a", "args_cutadapt"
    ]

    // Create a map from the input list 'it'
    def dit = [:]
    field_names.eachWithIndex { field_name, index ->
        if (index < it.size()) {
            dit[field_name] = it[index]
        } else {
            // Handle cases where optional fields might be missing at the end of the list
            // Or use default values from schema if necessary (more complex)
            dit[field_name] = null // Or some default
        }
    }

    def type = "SR"
    def r1_file = file(dit.path_r1, checkIfExists: true)
    def files = [r1_file]
    if (dit.path_r2 && !dit.path_r2.isEmpty()) {
        def r2_file = file(dit.path_r2, checkIfExists: true)
        files << r2_file
        type = "PE"
    }
    if (r1_file.toString().toLowerCase().endsWith("sfq")) {
        type = "sfq"
    }

    def meta = [
        "id"            : dit.sample_name,
        "reference"     : dit.reference,
        "strand"        : dit.strand,
        "is_3prime"     : dit.is_3prime as boolean,
        "frag_size"     : dit.frag_size as double,
        "frag_size_sd"  : dit.frag_size_sd as double,
        "read_type"     : type
    ]

    // Cutadapt specific meta construction
    def cutadapt_sample_args_list = []
    if (dit.r1_adapter && !dit.r1_adapter.isEmpty()) {
        cutadapt_sample_args_list.add("-a")
        cutadapt_sample_args_list.add("${dit.r1_adapter}")
    }
    if (type == "PE" && dit.r2_adapter && !dit.r2_adapter.isEmpty()) {
        cutadapt_sample_args_list.add("-A")
        cutadapt_sample_args_list.add("${dit.r2_adapter}")
    }
    
    boolean should_trim_poly_a = dit.containsKey('trim_poly_a') && dit.trim_poly_a != null ? (dit.trim_poly_a as boolean) : true
    if (should_trim_poly_a) {
        cutadapt_sample_args_list.add("--poly-a")
    }
  
    meta.args_cutadapt = cutadapt_sample_args_list.join(" ").trim() + (dit.args_cutadapt ?: "")


    return [meta, files]
}

def parse_reference(it) {
  def is_index = false
  def is_fasta = false
  def is_multi_fasta = false
  def ref_file
  def fasta_files = []
  
  if (it[2] && !it[2].isEmpty()) {
    ref_file = file(it[2], checkIfExists: true)
    is_index = true
  } else if (it[3] && !it[3].isEmpty()) {
    // Check if fasta field contains semicolons (multiple files)
    if (it[3].contains(';')) {
      is_multi_fasta = true
      is_fasta = true
      // Split by semicolon and check if each file exists
      def fasta_paths = it[3].split(';')
      fasta_paths.each { path ->
        def f = file(path.trim(), checkIfExists: true)
        fasta_files << f
      }
      ref_file = fasta_files
    } else {
      ref_file = file(it[3], checkIfExists: true)
      is_fasta = true
    }
  }
  
  // raise error if ref_file is not defined
  if (!ref_file) {
    throw new IllegalArgumentException("Reference file is not defined for ${it[0]} / ${it[1]} / ${it[2]} / ${it[3]}")
  }

  def meta = [
    "id": it[0], 
    "method": it[1],
    "is_index": is_index, 
    "is_fasta": is_fasta,
    "is_multi_fasta": is_multi_fasta
  ]
  return [meta, ref_file]
}

def make_quant_args(meta) {
  if (meta.method == "kallisto") {
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
    
    meta.kallisto_args = kallisto_args.join(" ")
  } else if (meta.method == "salmon") {
    def salmon_args = []

    if (meta.strand == "fr-stranded") {
      salmon_args << "--libType ISF"
    } else if (meta.strand == "rf-stranded") {
      salmon_args << "--libType ISR"
    } else {
      salmon_args << "--libType IU"
    }
    def has_frag_size = false
    if (meta.frag_size > 0) {
      salmon_args << "--fldMean " + meta.frag_size
      if (meta.frag_size_sd > 0) {
        salmon_args << "--fldSD " + meta.frag_size_sd
        has_frag_size = true
      }
    }

    if (meta.is_3prime) {
      salmon_args << "--noLengthCorrection"
      salmon_args << "--incompatPrior 0.0"
      if (!has_frag_size) {
        salmon_args << "--fldMean 250"
        salmon_args << "--fldSD 50"
      }
    }

    meta.args_salmon = salmon_args.join(" ")
  }
  
  return meta
}

include { PRIMARY } from './modules/subworkflows/primary/main.nf'
include { RNA_PREPROCESSING } from './modules/subworkflows/rna_preprocessing/main.nf'
include { KALLISTO_INDEX } from './modules/process/kallisto/index/main.nf'
include { SALMON_INDEX } from './modules/process/salmon/index/main.nf'
include { CONCAT_FILES } from './modules/process/concat_files/main.nf'
include { COMMIT_DIR } from './modules/process/commit_dir/main.nf'

workflow {

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
    log.info paramsHelp("nextflow run nexomis/rna-preprocessing --input </path/to/samplesheet> --reference </path/to/reference.csv> [args]")
    exit 0
  }
  validateParameters()
  log.info paramsSummaryLog(workflow)

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

    // Commit the kraken DB (directory or archive) to ensure it's a directory
    Channel.value(params.kraken2_db)
    | map { path ->
        def meta = [ id: "kraken_db" ] // Simple meta for the single DB item
        def input_item = file(path) // Let Nextflow stage the file/dir
        if (!input_item.exists()) {
            error "Provided kraken2_db path does not exist: ${path}"
        }
        return [ meta, input_item ]
      }
    | COMMIT_DIR
    | set { committedKrakenDb } // Output is tuple [meta, path_to_dir]

    taxDir = Channel.fromPath(params.tax_dir, type: 'dir', checkIfExists: true)
    numReads = Channel.value(params.num_reads_sample_qc)

    PRIMARY(rawInputs, committedKrakenDb, taxDir, numReads)

    trimmedInputs = PRIMARY.out.trimmed
  }
  // END PRIMARY

  // START RNA PREPROCESSING
  // Create a channel with reference methods
  Channel.fromSamplesheet("reference")
  | map {
      return parse_reference(it)
    }
  | set { reference }

  reference
  | map { [it[0].id, it[0].method] }
  | set { refMethods }

  // Join with trimmed inputs and add method to meta
  trimmedInputs
  | map { [it[0].reference, it] }
  | combine(refMethods, by: 0)
  | map { 
    def meta = it[1][0]
    meta.method = it[2]
    def meta_with_args = make_quant_args(meta)
    return [meta_with_args, it[1][1]]
  }
  | set { inputsForQuant }

  reference
  | branch {
    build: it[0].is_fasta
    existing: it[0].is_index
  }
  | set { referenceBranched }

  referenceBranched.build
  | branch {
    single: !it[0].is_multi_fasta
    multiple: it[0].is_multi_fasta
  }
  | set { refBuildBranched }

  // Concatenate multiple FASTA files if needed
  CONCAT_FILES(refBuildBranched.multiple)
  
  // Combine single and concatenated FASTA files
  singleFastaFiles = refBuildBranched.single
  concatenatedFastaFiles = CONCAT_FILES.out.fasta
  
  allFastaFiles = singleFastaFiles.concat(concatenatedFastaFiles)
  
  // Branch by method
  allFastaFiles
  | branch {
    kallisto: it[0].method == 'kallisto'
    salmon: it[0].method == 'salmon'
  }
  | set { refToBuild }

  KALLISTO_INDEX(refToBuild.kallisto)
  SALMON_INDEX(refToBuild.salmon)
  
  referenceIndexFinal = KALLISTO_INDEX.out.idx
    .concat(SALMON_INDEX.out.idx)
    .concat(referenceBranched.existing)

  RNA_PREPROCESSING(inputsForQuant, referenceIndexFinal)

  publish:
  PRIMARY.out.trimmed                         >> 'trimmed_and_filtered'
  PRIMARY.out.fastqc_trim_html                >> 'fastqc_trim'
  PRIMARY.out.fastqc_raw_html                 >> 'fastqc_raw'
  PRIMARY.out.multiqc_html                    >> 'multiqc'
  RNA_PREPROCESSING.out.kallisto_h5           >> 'kallisto'
  RNA_PREPROCESSING.out.kallisto_log          >> 'kallisto'
  RNA_PREPROCESSING.out.salmon_quant          >> 'salmon'
  RNA_PREPROCESSING.out.multiqc               >> 'align_multiqc'
}

output {
  'trimmed_and_filtered' {
    enabled params.save_trimmed_reads
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
  'salmon' {
    path 'salmon'
  }
  'align_multiqc' {
    path 'qc/multiqc/align'
  }
}
