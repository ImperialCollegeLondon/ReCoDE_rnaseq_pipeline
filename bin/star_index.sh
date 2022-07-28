#!/bin/bash
# 
# This script indexes an input genome using STAR
#
# The script takes the following inputs:
#
# $1 = directory to save indices (e.g. results/1_basic_pipeline/e_star_index/)
# $2 = genome fasta files (e.g. data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna)
# $3 = genome annotation file file (e.g. data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf)
# $4 = scale, must be reduced for smaller genomes
# $5 = sparsity (can be set to non-zero values to reduce the memory STAR needs)
# $6 = number of threads for STAR to use

# genome files must be unzipped for STAR indexing
gzip -cfdk "$2" > "$1/genome.fna"
gzip -cfdk "$3" > "$1/annotation.gtf"

# generate genome indices
# can use --genomeSAsparseD as an argument if STAR uses too much memory
STAR --runMode genomeGenerate \
     --genomeDir "$1" \
     --genomeFastaFiles "$1/genome.fna" \
     --sjdbGTFfile "$1/annotation.gtf" \
     --genomeSAindexNbases "$4" \
     --genomeSAsparseD "$5" \
     --runThreadN "$6"
