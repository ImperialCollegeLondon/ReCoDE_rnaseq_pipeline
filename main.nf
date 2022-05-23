#!/usr/bin/env nextflow

// enable module feature
nextflow.enable.dsl = 2

include { PROCESS_RNASEQ } from './workflows/nextflow_pipeline'

// run workflow in workflows/
workflow {
    PROCESS_RNASEQ()
}
