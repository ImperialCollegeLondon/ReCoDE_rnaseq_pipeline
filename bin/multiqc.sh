#!/bin/bash
# 
# This script runs multiqc on up to three directories
#
# The script takes the following inputs:
#
# $1 = the results directory (e.g. results/1_basic_pipeline/b_multiqc)
# $2 = directory 1 (e.g. results/1_basic_pipeline/a_fastqc)
# $3 = directory 2 (e.g. results/1_basic_pipeline/f_align)
# $4 = directory 3 (e.g. results/1_basic_pipeline/g_count)

# group together fastqc reports
multiqc -f -o "$1" "$2" "$3" "$4"
