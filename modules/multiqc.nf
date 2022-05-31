process MULTIQC {

    label "short"
    publishDir "$params.outdir/multiqc/", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::multiqc=1.12 conda-forge::spectra=0.0.11" : null)
    container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"

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
