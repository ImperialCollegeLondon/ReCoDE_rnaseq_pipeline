#!/bin/bash

# $1 = the results directory (e.g. results/1_basic_pipeline/a_fastqc)
# $2 = the path of the fastq file (e.g. data/fastq/SRR391535.fastq.gz)
# $3 = temporary directory (e.g. results/1_basic_pipeline/a_fastqc/SRR391535)

# make the temporary directory
mkdir $3

# run fastqc on the specified fastq file
fastqc -o $1 -d $3 $2
