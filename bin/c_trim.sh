
# $1 = the results directory (e.g. results/1_basic_pipeline/c_trim)
# $2 = the path of the fastq file (e.g. data/fastq/SRR391535.fastq.gz)

# run fastqc on the specified fastq file
trim_galore --fastqc -o $1 $2
# --fastqc_args "<ARGS>"
