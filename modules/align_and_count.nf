process ALIGN_AND_COUNT {

    label "star_align"
    publishDir "$params.outdir/f_align_and_count/", mode: params.publish_dir_mode

    input:
    tuple val(accession), path(fastq)
    path indexed_genome
    path annotation

    output:
    tuple val(accession), path("*.counts"), emit: counts
    tuple val(accession), path("*Log.final.out"), emit: log_final

    script:
    """
    $baseDir/bin/align_and_count.sh \
        "$indexed_genome" \
        "$fastq" \
        "$accession" \
        "$task.cpus" \
        "$annotation"
    """
}
