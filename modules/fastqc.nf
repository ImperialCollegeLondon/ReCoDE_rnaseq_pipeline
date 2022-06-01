process FASTQC {

    label "long"
    publishDir "$params.outdir/a_fastqc/", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)

    // use a custom container rather than biocontainers
    // container "quay.io/biocontainers/fastqc:0.11.9--0"
    container "jackgisby/fastqc"

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
