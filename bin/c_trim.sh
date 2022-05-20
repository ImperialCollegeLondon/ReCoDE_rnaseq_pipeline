
# $1 = the results directory (e.g. results/1_basic_pipeline/c_trim)
# $2 = the path of the fastq file (e.g. data/fastq/SRR391535.fastq.gz)
# $3 = sample name (e.g. SRR391535)

# make temporary directory
mkdir $1/$3

# run fastqc on the specified fastq file  
trim_galore --fastqc --fastqc_args "-d $1/$3" -o $1 $2

rmdir $1/$3
