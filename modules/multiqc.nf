process MULTIQC {

    label "short"
    publishDir "$params.outdir/multiqc/", mode: params.publish_dir_mode

    conda "bioconda::multiqc=1.12 conda-forge::spectra=0.0.11"

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
