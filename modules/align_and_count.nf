process ALIGN_AND_COUNT {

    label "star_align"

    input:
    tuple val(accession), path(fastq)
    path indexed_genome
    path unzipped_annotation

    output:
    tuple val(accession), path("*.counts"), emit: counts
    path "*Log.final.out",                  emit: log_final

    script:
    """
    $baseDir/bin/align_and_count.sh \
        $indexed_genome \
        $fastq \
        $accession \
        $task.cpus \
        $unzipped_annotation
    """
}
