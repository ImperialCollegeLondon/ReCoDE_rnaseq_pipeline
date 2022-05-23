// import local modules
include { GET_INPUT_VARIANTS } from "$baseDir/modules/get_input_variants.nf"

// main pipeline
workflow PROCESS_RNASEQ {

    // get input files
    input_fastq = Channel.fromPath("${params.fastqc_dir}/*.fastq.gz")

    // apply QC to the raw reads
    FASTQC(input_fastq)
    MULTIQC()

    // // trim the raw reads and re-apply QC
    // TRIM()
    // MULTIQC()

    // // index the genome for alignment
    // STAR_INDEX()

    // // align reads and get gene counts
    // ALIGN_AND_COUNT()
    // MULTIQC()
}
