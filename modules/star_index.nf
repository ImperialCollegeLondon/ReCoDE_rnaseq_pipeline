process STAR_INDEX {

    label "star_index"
    publishDir "$params.outdir/e_star_index/", mode: params.publish_dir_mode

    input:
    path genome
    path annotation

    output:
    path "*", emit: indexed_genome
    path "*.gtf", emit: unzipped_annotation

    script:
    """
    $baseDir/bin/star_index.sh \
        "./" \
        $genome \
        $annotation \
        $params.sparse_d \
        $task.cpus
    """
}
