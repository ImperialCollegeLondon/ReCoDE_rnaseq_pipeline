
# $1 = directory to save indices (e.g. results/1_basic_pipeline/e_star_index/)
# $2 = genome fasta files (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.fna)
# $3 = genome annotation file file (e.g. results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf)
# $4 = sparsity (can be set to non-zero values to reduce the memory STAR needs)
# $5 = number of threads for STAR to use

# genome files must be unzipped for STAR indexing
gzip -cfdk $2 > $1/GCF_000004515.6_Glycine_max_v4.0_genomic.fna
gzip -cfdk $3 > $1/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

# generate genome indices
# can use --genomeSAsparseD as an argument if STAR uses too much memory
STAR --runMode genomeGenerate \
     --genomeDir $1 \
     --genomeFastaFiles $1/GCF_000004515.6_Glycine_max_v4.0_genomic.fna \
     --sjdbGTFfile $1/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf \
     --genomeSAindexNbases 13 \
     --genomeSAsparseD $4 \
     --runThreadN $5

# these lines needed to be removed for htseq to run
sed -i "1326025d;1326521d" $1/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf
