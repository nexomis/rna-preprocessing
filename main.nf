#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

log.info """
    |            #################################################
    |            #    _  _                             _         #
    |            #   | \\| |  ___  __ __  ___   _ __   (_)  __    #
    |            #   | .` | / -_) \\ \\ / / _ \\ | '  \\  | | (_-<   #
    |            #   |_|\\_| \\___| /_\\_\\ \\___/ |_|_|_| |_| /__/   #
    |            #                                               #
    |            #################################################
    |
    | DE-Preprocessing: Estimate the expression level of features (genes or transcripts) and perform the associated QC. Starting from cleaned reads and one pre-built reference database in standard format (cf. Readme.md).
    | 
    |""".stripMargin()

if (params.help) {
  log.info paramsHelp("nextflow run nexomis/DE-Preprocessing --input_dir </path/to/dir/contains/fastq> --kallisto_idx </path/to/kallisto/txome/file.index> [args]")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

include { RNA_ALIGN_KALLISTO } from './modules/subworkflows/rna_align_kallisto/main.nf'

workflow {

  // formatted reads channels
  inputDir = Channel.fromPath(params.input_dir, type:'dir', checkIfExists: true)
  kallistoIdx = Channel.fromPath(params.kallisto_idx, type:'file', checkIfExists: true) | first()  // check convention: in tuple with 'name' at first position ?
  readsOrientation = Channel.value(params.reads_orientation)

  // execution
  RNA_ALIGN_KALLISTO(inputDir, kallistoIdx, readsOrientation)
}
