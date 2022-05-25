// import local modules
include { FASTQC } from "$baseDir/modules/fastqc.nf"
include { MULTIQC } from "$baseDir/modules/multiqc.nf"
include { TRIM } from "$baseDir/modules/trim.nf"
include { STAR_INDEX } from "$baseDir/modules/star_index.nf"
include { ALIGN_AND_COUNT } from "$baseDir/modules/align_and_count.nf"

// main pipeline
workflow PROCESS_RNASEQ {

    // get input fastq and their accession IDs
    input_fastq = Channel.fromPath("${params.fastqc_dir}/*.fastq.gz").map {
        tuple( it.name.split('\\.')[0], it )
    }

    // apply QC to the raw reads
    FASTQC(input_fastq)

    // option to only do QC, skipping trimming and alignment
    if (params.trim) {

        // trim the raw reads and re-apply QC
        TRIM(input_fastq)

        // define channels
        ch_trimmed_fastqc_html = TRIM.out.trimmed_fastqc_html
        ch_trimmed_fastqc_zip  = TRIM.out.trimmed_fastqc_zip
        ch_trimming_report     = TRIM.out.trimming_report
        ch_fastq_to_align      = TRIM.out.trimmed_fastq

    } else {

        // if we don't trim, we must align the raw reads
        ch_fastq_to_align = input_fastq

        // trimming skipped, define channels as empty for input to multiqc
        ch_trimmed_fastqc_html = Channel.empty()
        ch_trimmed_fastqc_zip  = Channel.empty()
        ch_trimming_report     = Channel.empty()
    }

    // option to only trim and QC, skipping alignment
    if (params.align) {

        // index the genome for alignment
        STAR_INDEX(file("$params.genome_fasta"), file("$params.genome_gtf"))

        // align reads and get gene counts
        ALIGN_AND_COUNT(
            ch_fastq_to_align, 
            STAR_INDEX.out.indexed_genome, 
            STAR_INDEX.out.annotation
        )

        // define channels for multiqc
        ch_star_log = ALIGN_AND_COUNT.out.log_final
        ch_counts = ALIGN_AND_COUNT.out.counts

    } else {

        // define empty channels for multiqc
        ch_star_log = Channel.empty()
        ch_counts = Channel.empty()
    }

    MULTIQC(
        FASTQC.out.fastqc_zip.collect().ifEmpty([]),
        ch_trimmed_fastqc_zip.collect().ifEmpty([]),
        ch_trimming_report.collect().ifEmpty([]),
        ch_star_log.collect().ifEmpty([]),
        ch_counts.collect().ifEmpty([])
    )
}
