
# the variable $1 is used to get the argument we passed to this script
echo $1

# we get the value of the row specified by $1 here
fastq_srr="$(head -"$1" data/files.txt | tail -1)"
echo $fastq_srr

# so, this is the location of the .fastq file we need
fastq_location="data/$fastq_srr.fastq.gz"
echo $fastq_location

# run fastqc on the specified fastq file
fastqc $fastq_location -o results/1_basic_pipeline/a_fastqc
