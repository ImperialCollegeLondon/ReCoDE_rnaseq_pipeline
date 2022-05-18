
# stop running the script if there are errors
set -e

# download the data (skip this step if already downloaded)
if [ ! -f data/files.txt ]; then
    data/get_data.sh
fi

# get number of samples to run
readarray -t SAMPLE_SRR < data/files.txt

########################################## run the pipeline steps

# make the results directory if needed
if [ ! -e  results/1_basic_pipeline/a_fastqc ]; then
    mkdir results/1_basic_pipeline/a_fastqc
fi

# for each sample
for s in "${SAMPLE_SRR[@]}"; do

    # run fastqc on raw fastq
    scripts/1_basic_pipeline/a_fastqc.sh \
        results/1_basic_pipeline/a_fastqc \  # $1 = the results directory
        data/fastq/$s.fastq.gz               # $2 = the path of the fastq file
done

# combine the fastqc results
scripts/1_basic_pipeline/b_multiqc.sh \
        results/1_basic_pipeline/b_multiqc \  # $1 = the results directory
        results/1_basic_pipeline/a_fastqc     # $2 = directory of the fastqc output

for s in "${SAMPLE_SRR[@]}"; do

    # trim the fastq files
    scripts/1_basic_pipeline/c_trim.sh \
        results/1_basic_pipeline/c_trim \  # $1 = the results directory
        data/fastq/$s.fastq.gz             # $2 = the path of the fastq file
done

# combine the fastqc results generated for the trimmed fastq files
scripts/1_basic_pipeline/d_trim_multiqc.sh \
        results/1_basic_pipeline/d_trim_multiqc \  # $1 = the results directory 
        results/1_basic_pipeline/c_trim            # $2 = directory of the fastqc output

# genome files must be unzipped for STAR indexing
gunzip -fdk data/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz -C results/1_basic_pipeline/e_star_index/
gunzip -fdk data/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz -C results/1_basic_pipeline/e_star_index/

# index the genome using STAR
scripts/1_basic_pipeline/e_star_index.sh \
        results/1_basic_pipeline/e_star_index/ \                                              # $1 = directory to save indices
        results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.fna \  # $2 = genome fasta files
        results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf    # $3 = genome annotation file file

# remove unzipped fasta
rm results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.fna

for s in "${SAMPLE_SRR[@]}"; do

    # perform alignment using STAR, providing the directory of the indexed genome
    scripts/1_basic_pipeline/f_star_align.sh \
        results/1_basic_pipeline/e_star_index \        # $1 = the location of the indexed genome
        results/1_basic_pipeline/c_trim/$s.fastq.gz \  # $2 = the location of the trimmed fastq file
        results/1_basic_pipeline/f_star_align/$s       # $3 = the prefix of the output file

    # generate a counts matrix from the alignment
    scripts/1_basic_pipeline/g_htseq_count.sh \
        results/1_basic_pipeline/f_star_align/$s.bam \                                        # $1 = the alignment file
        results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf \  # $2 = the annotation file
        results/1_basic_pipeline/g_htseq_count/$s                                             # $3 = the prefix of the output file
done

# remove unzipped gtf
rm results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

# combine the counts for each sample into a single matrix
scripts/1_basic_pipeline/h_combine_counts.sh \
    results/1_basic_pipeline/g_htseq_count/ \                      # $1 = directory containing counts
    results/1_basic_pipeline/h_combine_counts/combined_counts.txt  # $2 = output file path
