process TRIM {

    label "long"
    publishDir "$params.outdir/c_trim/", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::trim-galore=0.6.7" : null)
    container "quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"

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

    $baseDir/bin/trim.sh \
        "./" \
        "$fastq" \
        "./trim_tmp"
    """
}
