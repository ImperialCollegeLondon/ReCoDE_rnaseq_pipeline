process ALIGN {

    label "star_align"
    publishDir "$params.outdir/f_align/", mode: params.publish_dir_mode

    conda "bioconda::star=2.7.10a"

    input:
    tuple val(accession), path(fastq)
    path indexed_genome

    output:
    tuple val(accession), path("*Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    tuple val(accession), path("*Log.final.out"), emit: log_final

    script:
    """
    $baseDir/bin/align.sh \
        "$indexed_genome" \
        "$fastq" \
        "$accession" \
        "$task.cpus"
    """
}

process COUNT {

    label "long"
    publishDir "$params.outdir/g_count/", mode: params.publish_dir_mode

    conda "bioconda::htseq=0.11.3"

    input:
    tuple val(accession), path(bam)
    path annotation

    output:
    tuple val(accession), path("*.counts"), emit: counts

    script:
    """
    $baseDir/bin/align_and_count.sh \
        "$bam" \
        "$annotation" \
        "$accession"
    """
}
