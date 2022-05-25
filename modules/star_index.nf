process STAR_INDEX {

    label "star_index"
    publishDir "$params.outdir/e_star_index/", mode: params.publish_dir_mode

    input:
    path genome
    path annotation

    output:
    path "indexed_genome", emit: indexed_genome
    path "indexed_genome/SAindex", emit: sa_index
    path "*.gtf", emit: annotation

    script:
    """
    mkdir indexed_genome

    $baseDir/bin/star_index.sh \
        "indexed_genome" \
        "$genome" \
        "$annotation" \
        "$params.sparse_d" \
        "$task.cpus"

    mv indexed_genome/*gtf ./
    """
}
