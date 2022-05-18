
# $1 = directory to save indices (e.g. results/1_basic_pipeline/e_star_index/)
# $2 = genome fasta files (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.fna)
# $3 = genome annotation file file (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf)

# generate genome indices
# can use --genomeSAsparseD 2 as an argument if STAR uses too much memory
STAR --runMode genomeGenerate \
     --genomeDir $1 \
     --genomeFastaFiles $2 \
     --sjdbGTFfile $3
