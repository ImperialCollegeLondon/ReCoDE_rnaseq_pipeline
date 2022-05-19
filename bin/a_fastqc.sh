
# $1 = the results directory (e.g. results/1_basic_pipeline/a_fastqc)
# $2 = the path of the fastq file (e.g. data/fastq/SRR391535.fastq.gz)

# run fastqc on the specified fastq file
fastqc -o $1 $2
