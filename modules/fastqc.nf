process FASTQC {

    label "long"
    publishDir "$params.outdir/a_fastqc/", mode: params.publish_dir_mode

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
