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
  log.info paramsHelp("nextflow run nexomis/DE-Preprocessing [args]")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

include { PARSE_SEQ_DIR_UNSPRING } from './modules/subworkflows/parse_seq_dir_unspring/main.nf'
include { KALLISTO_QUANT } from './modules/process/kallisto/quant/main.nf'

workflow {
  // formatted reads channels
  inputDir = Channel.fromPath(params.input_dir, type:'dir', checkIfExists: true)
  reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

  // kallisto quant
  kallisto_idx = Channel.fromPath(params.kallisto_idx, type:'file', checkIfExists: true)
  KALLISTO_QUANT(reads, kallisto_idx.first())
}