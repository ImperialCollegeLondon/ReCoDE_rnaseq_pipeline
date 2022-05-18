
# $1 = the alignment file (e.g. results/1_basic_pipeline/f_star_align/SRR391535.bam)
# $2 = the annotation file (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf)
# $3 = the prefix of the output file (e.g. results/1_basic_pipeline/g_htseq_count/SRR391535)

htseq-count -s no -r pos -f bam $1 $2 > $3.counts
