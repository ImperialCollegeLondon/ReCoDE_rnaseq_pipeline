
# $1 = the location of the indexed genome (e.g. results/1_basic_pipeline/e_star_index)
# $2 = the location of the trimmed fastq file (e.g. results/1_basic_pipeline/c_trim/SRR391535_trimmed.fq.gz)
# $3 = the prefix of the output file (e.g. results/1_basic_pipeline/f_star_align/SRR391535)
# $4 = the number of threads for STAR to use
# $5 = the annotation file (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf) 

# unzip fastq for alignment
gzip -cfdk $2 > $3.fastq

# perform alignment
STAR --genomeDir $1 \
     --readFilesIn $3.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $3 \
     --runThreadN $4

# remove unzipped fastq
rm $3.fastq

# index the alignment for htseq
samtools index $3Aligned.sortedByCoord.out.bam

# use alignment to get counts
htseq-count -s no -r pos -f bam $3Aligned.sortedByCoord.out.bam $5 > $3.counts
