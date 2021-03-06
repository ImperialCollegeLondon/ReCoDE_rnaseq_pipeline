#!/bin/bash
# 
# This script quantifies the reads from the alignment file
#
# The script takes the following inputs:
#
# $1 = the input bam file (e.g. results/1_basic_pipeline/f_align/SRR391535Aligned.sortedByCoord.out.bam)
# $2 = the annotation file (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf) 
# $3 = the prefix of the output file (e.g. results/1_basic_pipeline/g_count/SRR391535)

# use alignment to get counts
htseq-count -s no -r pos -f bam "$1" "$2" > "$3.counts"
