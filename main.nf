#!/usr/bin/env nextflow

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
  def meta = ["id": it[0], "strand": it[3], "is_3prime": it[4],
      "frag_size": it[5],"frag_size_sd": it[6],
      "force_kallisto_single_overhang": it[7]
    ]
  def r1_file = file(it[1])
  if (it[2] == "") {
    if (r1_file.getExtension() == "spring") {
      meta.read_type = "spring"
    } else {
      meta.read_type = "SR"
    }
    return [meta, [r1_file] ]
  } else {
    def r2_file = file(it[2])
    meta.read_type = "PE"
    return [meta, [r1_file, r2_file] ]
  }
}

def make_kallisto_args(meta) {
  kallisto_args = []
  if (meta.strand == "fr-stranded") {
    kallisto_args << "--fr-stranded"
  } else if (meta.strand == "rf-stranded") {
    kallisto_args << "--rf-stranded"
  }
  if (meta.is_3prime || meta.force_kallisto_single_overhang) {
    kallisto_args << "--single-overhang"
  }
  if (meta.frag_size > 0) {
    kallisto_args << "--fragment-length " + meta.frag_size
  } else {
    if ( meta.read_type == "SR" ) {
      kallisto_args << "--fragment-length 250"
    }
  }
  if (meta.frag_size_sd > 0) {
    kallisto_args << "--fragment-length " + meta.frag_size_sd
  } else {
    if ( meta.read_type == "SR" ) {
      kallisto_args << "--sd 50"
    }
  }

  return kallisto_args.join(" ")
}

process ALIGN_MULTIQC {
  container 'multiqc/multiqc:v1.21' // version above are bugged

  label 'cpu_x1'
  label 'mem_8G'

  input:
  path kallisto_logs, stageAs: 'kallisto_logs/*', arity: '1..*'
  path conf_yml, arity: 1, stageAs: "multiqc_config.yaml"

  output:
  path("align_multiqc.html", arity: 1)

  script:
  """
  #!/usr/bin/bash

  multiqc -c $conf_yml --no-data-dir -n align_multiqc.html .

  """

  stub:
  """
  #!/usr/bin/bash

  touch primary_multiqc.html

  """
}

include {SPRING_DECOMPRESS} from './modules/process/spring/decompress/main.nf'
include {KALLISTO_QUANT} from './modules/process/kallisto/quant/main.nf'
include {PRIMARY_FROM_READS} from './modules/subworkflows/primary/from_reads/main.nf'

workflow {

  // START PARSING SAMPLE SHEET
  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | branch {
    spring: it[0].read_type == "spring"
    fastq: it[0].read_type != "spring"
  }
  | set { inputs }

  SPRING_DECOMPRESS(inputs.spring)
  | map {
    if (it[1].size()==1) {
      it[0].read_type = "SR"
    } else {
      it[0].read_type = "PE"
    }
    return it
  }
  | concat(inputs.fastq)
  | set {rawInputs}
  // END PARSING SAMPLE SHEET

  // START PRIMARY
  if (params.skip_primary) {
    trimmedInputs = rawInputs
  } else {
    if ( params.kraken2_db == null ) {
      error "kraken2_db argument required for primary analysis"
    }
    Channel.fromPath(params.kraken2_db, type: "dir", checkIfExists: true)
    | collect
    | set {dbPathKraken2}
    
    PRIMARY_FROM_READS(rawInputs, dbPathKraken2)
    
    trimmedInputs = PRIMARY_FROM_READS.out.trimmed
  }
  // END PRIMARY

  // START KALLISTO
  trimmedInputs
  | map {
    it[0].kallisto_args = make_kallisto_args(it[0])
    return it
  }
  | set {inputsForKallisto}

  kallistoIndex = Channel.fromPath(params.kallisto_idx).collect()
  
  KALLISTO_QUANT(inputsForKallisto, kallistoIndex)

  multiqcYml = Channel.fromPath(projectDir + "/files/align_multiqc.yml")

  KALLISTO_QUANT.out.log
  | map {it[1]}
  | collect
  | set {kallistoLogs}
  // END KALLISTO

  ALIGN_MULTIQC(
    kallistoLogs,
    multiqcYml
  )

}
