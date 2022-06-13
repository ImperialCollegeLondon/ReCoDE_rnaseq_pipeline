process ALIGN {

    label "star_align"
    publishDir "$params.outdir/f_align/", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::star=2.7.10a bioconda::samtools=1.6" : null)
    container "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0"

    input:
    tuple val(accession), path(fastq)
    path indexed_genome

    output:
    tuple val(accession), path("*Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    path "*Log.final.out", emit: log_final

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

    conda (params.enable_conda ? "bioconda::htseq=0.11.3" : null)
    container "biocontainers/htseq:v0.11.2-1-deb-py3_cv1"
    
    input:
    tuple val(accession), path(bam)
    path annotation

    output:
    path "*.counts", emit: counts

    script:
    """
    $baseDir/bin/count.sh \
        "$bam" \
        "$annotation" \
        "$accession"
    """
}
