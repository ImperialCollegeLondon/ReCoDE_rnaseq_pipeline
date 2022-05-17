
# get the location of the fastq file
fastq_srr="$(head -"$1" data/files.txt | tail -1)"
fastq_location="data/$fastq_srr.fastq.gz"

# run fastqc on the specified fastq file
trim_galore --fastqc -o results/1_basic_pipeline/b_trim $fastq_location 

# --cores INT
