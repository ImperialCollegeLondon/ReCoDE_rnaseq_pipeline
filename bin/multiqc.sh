
# $1 = the results directory (e.g. results/1_basic_pipeline/b_multiqc)
# $2 = directory 1 (e.g. results/1_basic_pipeline/a_fastqc)
# $3 = directory 2 (e.g. results/1_basic_pipeline/f_align_and_count)

# group together fastqc reports
multiqc -f -o $1 $2 $3
