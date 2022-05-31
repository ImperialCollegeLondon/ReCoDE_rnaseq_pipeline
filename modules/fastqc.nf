process FASTQC {

    label "long"
    publishDir "$params.outdir/a_fastqc/", mode: params.publish_dir_mode

    // setup conda/containers, based on nf-core/rnaseq
    // if conda is enabled, use package from bioconda
    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)

    // get docker/singularity image if docker/singularity is being used
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

    input:
    tuple val(accession), path(fastq)

    output:
    path "*_fastqc.html", emit: fastqc_html
    path "*_fastqc.zip",  emit: fastqc_zip

    script:
    """
    mkdir fastqc_tmp

    $baseDir/bin/fastqc.sh \
        "./" \
        "$fastq" \
        "./fastqc_tmp"
    """
}
