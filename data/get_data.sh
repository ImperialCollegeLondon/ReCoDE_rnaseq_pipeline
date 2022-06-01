#!/bin/bash
#
# This script downloads the full reads and genome to be analysed by the pipeline
# For a smaller scale test, the repository comes with a subset of these data

# a vector of the sequence read archive identifiers for the RNA-seq samples
readarray -t SAMPLE_SRR < data/files.txt

# load fastq-dump command on the cluster
# module load sra-toolkit/2.8.1

# create a file to list the IDs and file locations
rm -f data/files.txt
touch data/files.txt

# for each sample
for srr in ${SAMPLE_SRR[@]}; do

  # get the sample's file from the sequence read archive
  fastq-dump --gzip -O data/fastq/ "$srr"

done

# get the genome and annotation files from NCBI
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000004515.6

mkdir data/genome

wget -O data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz

wget -O data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz
