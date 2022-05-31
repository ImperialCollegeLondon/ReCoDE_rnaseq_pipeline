process MULTIQC {

    label "short"
    publishDir "$params.outdir/multiqc/", mode: params.publish_dir_mode

    // setup conda/containers, based on nf-core/rnaseq
    // if conda is enabled, use package from bioconda
    conda (params.enable_conda ? "bioconda::multiqc=1.12 conda-forge::spectra=0.0.11" : null)

    // get docker/singularity image if docker/singularity is being used
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0' :
            'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path ("fastqc/*")
    path ("trimmed/fastqc/*")
    path ("trimmed/*")
    path ("star/*")
    path ("htseq/*")

    output:
    path "*multiqc_report.html"
    path "*_data"

    script:
    """
    $baseDir/bin/multiqc.sh \
        "./" \
        "--config $multiqc_config" \
        "./"
    """
}
