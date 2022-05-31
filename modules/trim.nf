process TRIM {

    label "long"
    publishDir "$params.outdir/c_trim/", mode: params.publish_dir_mode

    // setup conda/containers, based on nf-core/rnaseq
    // if conda is enabled, use package from bioconda
    conda (params.enable_conda ? "bioconda::trim-galore=0.6.7" : null)

    // get docker/singularity image if docker/singularity is being used
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    input:
    tuple val(accession), path(fastq)

    output:
    tuple val(accession), path("*_trimmed.fq.gz"), emit: trimmed_fastq
    path "*_trimmed_fastqc.html",                  emit: trimmed_fastqc_html
    path "*_trimmed_fastqc.zip",                   emit: trimmed_fastqc_zip
    path "*.fastq.gz_trimming_report.txt",         emit: trimming_report

    script:
    """
    mkdir trim_tmp

    conda list
    trim_galore -v

    $baseDir/bin/trim.sh \
        "./" \
        "$fastq" \
        "./trim_tmp"
    """
}
