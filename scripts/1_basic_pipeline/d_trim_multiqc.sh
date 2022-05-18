
# $1 = the results directory (e.g. results/1_basic_pipeline/d_trim_multiqc)
# $2 = directory of the fastqc output (e.g. results/1_basic_pipeline/c_trim)

# group together fastqc reports
multiqc -o $1 $2
