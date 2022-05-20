
# $1 = the results directory (e.g. results/1_basic_pipeline/g_final_multiqc)
# $2 = directory of the fastqc output (e.g. results/1_basic_pipeline/c_trim)
# $3 = directory of the alignments and counts (e.g. results/1_basic_pipeline/f_align_and_count)

# group together fastqc reports with the alignment and counts
multiqc -o $1 $2 $3
