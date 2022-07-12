#!/bin/bash
#
# This script downloads the full reads and genome to be analysed by the pipeline
# For a smaller scale test, the repository comes with a subset of these data

# a vector of the sequence read archive identifiers for the RNA-seq samples
readarray -t SAMPLE_SRR < data/files.txt

# load fastq-dump command on the cluster
# module load sra-toolkit/2.8.1

# for each sample
# for srr in ${SAMPLE_SRR[@]}; do

#   # get the sample's file from the sequence read archive
#   fastq-dump --gzip -O data/fastq/ "$srr"

# done

# get the genome and annotation files from NCBI
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000004515.6

mkdir data/genome
SOYBEAN_GENOME_ACCESSION=GCF_000004515.6_Glycine_max_v4.0_genomic

wget -O "data/genome/${SOYBEAN_GENOME_ACCESSION}.fna.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/${SOYBEAN_GENOME_ACCESSION}.fna.gz"

wget -O "data/genome/${SOYBEAN_GENOME_ACCESSION}.gtf.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/${SOYBEAN_GENOME_ACCESSION}.gtf.gz"

# there are a couple of lines in the soybean genome annotation that are misformatted
# the following code removes these lines so that there aren't any issues later

# decompress
gzip -cfd "data/genome/${SOYBEAN_GENOME_ACCESSION}.gtf.gz" > "data/genome/${SOYBEAN_GENOME_ACCESSION}.gtf"

# remove lines
sed -i "1326025d;1326521d" "data/genome/${SOYBEAN_GENOME_ACCESSION}.gtf"

# recompress
gzip -cf "data/genome/${SOYBEAN_GENOME_ACCESSION}.gtf" > "data/genome/${SOYBEAN_GENOME_ACCESSION}.gtf.gz"
