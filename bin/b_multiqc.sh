
# $1 = the results directory (e.g. results/1_basic_pipeline/b_multiqc)
# $2 = directory of the fastqc output (e.g. results/1_basic_pipeline/a_fastqc)

# group together fastqc reports
multiqc -o $1 $2
