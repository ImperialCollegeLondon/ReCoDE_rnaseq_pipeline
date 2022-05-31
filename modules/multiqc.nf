process MULTIQC {

    label "short"
    publishDir "$params.outdir/multiqc/", mode: params.publish_dir_mode

    conda "bioconda::multiqc=1.12"

    input:
    path ("fastqc/*")
    path ("trimmed/fastqc/*")
    path ("trimmed/*")
    path ("star/*")
    path ("htseq/*")

    output:
    path "multiqc_output/*multiqc_report.html"
    path "multiqc_output/*_data"

    script:
    """
    ls
    
    $baseDir/bin/multiqc.sh \
        "./multiqc_output" \
        "./"
    """
}
