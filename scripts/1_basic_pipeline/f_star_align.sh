
# $1 = the location of the indexed genome (e.g. results/1_basic_pipeline/e_star_index)
# $2 = the location of the trimmed fastq file (e.g. results/1_basic_pipeline/c_trim/SRR391535.fastq.gz)
# $3 = the prefix of the output file (e.g. results/1_basic_pipeline/f_star_align/SRR391535)

mkdir results/1_basic_pipeline/f_star_align

# perform alignment
STAR --genomeDir $1 \
     --readFilesCommand gunzip -c \
     --readFilesIn $2 \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $3
