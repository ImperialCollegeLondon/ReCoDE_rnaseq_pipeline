process STAR_INDEX {

    label "star_index"
    publishDir "$params.outdir/e_star_index/", mode: params.publish_dir_mode

    // setup conda/containers, based on nf-core/rnaseq
    // if conda is enabled, use package from bioconda
    conda (params.enable_conda ? "bioconda::star=2.7.10a" : null)

    // get docker/singularity image if docker/singularity is being used
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0' }"

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
