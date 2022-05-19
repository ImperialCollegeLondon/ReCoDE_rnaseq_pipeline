
# stop running the script if there are errors
set -e

# download the data (skip this step if already downloaded)
if [ ! -f data/files.txt ]; then
    data/get_data.sh
fi

# get number of samples to run
readarray -t SAMPLE_SRR < data/files.txt

########################################## run the pipeline steps

# optional: remove results from previous run
# rm -rf results/1_basic_pipeline
# mkdir results/1_basic_pipeline

# function creates a results folder if it doesn't already exist
create_folder () {
    if [ ! -e  results/1_basic_pipeline/$1 ]; then
        mkdir results/1_basic_pipeline/$1
    fi
}

create_folder "a_fastqc"

# for each sample
for s in "${SAMPLE_SRR[@]}"; do

    # run fastqc on raw fastq
    scripts/1_basic_pipeline/a_fastqc.sh \
        results/1_basic_pipeline/a_fastqc \
        data/fastq/$s.fastq.gz
done

# combine the fastqc results
scripts/1_basic_pipeline/b_multiqc.sh \
        results/1_basic_pipeline/b_multiqc \
        results/1_basic_pipeline/a_fastqc

for s in "${SAMPLE_SRR[@]}"; do

    # trim the fastq files
    scripts/1_basic_pipeline/c_trim.sh \
        results/1_basic_pipeline/c_trim \
        data/fastq/$s.fastq.gz
done

# combine the fastqc results generated for the trimmed fastq files
scripts/1_basic_pipeline/d_trim_multiqc.sh \
        results/1_basic_pipeline/d_trim_multiqc \
        results/1_basic_pipeline/c_trim

create_folder "e_star_index"

# index the genome using STAR
scripts/1_basic_pipeline/e_star_index.sh \
        results/1_basic_pipeline/e_star_index/ \
        data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz \
        data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz \
        3 \
        6

create_folder "f_align_and_count"

for s in "${SAMPLE_SRR[@]}"; do


    # perform alignment using STAR, providing the directory of the indexed genome
    scripts/1_basic_pipeline/f_align_and_count.sh \
        results/1_basic_pipeline/e_star_index \
        results/1_basic_pipeline/c_trim/"$s"_trimmed.fq.gz \
        results/1_basic_pipeline/f_align_and_count/$s \
        6 \
        results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

    # generate a counts matrix from the alignment
    scripts/1_basic_pipeline/g_htseq_count.sh \
        results/1_basic_pipeline/f_align_and_count/"$s"Aligned.sortedByCoord.out.bam \
        
        results/1_basic_pipeline/f_align_and_count/$s
done

# remove unzipped gtf
rm results/1_basic_pipeline/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

# use multiqc to assess the alignment and counts
scripts/1_basic_pipeline/g_final_multiqc.sh \
    results/1_basic_pipeline/g_final_multiqc \
    results/1_basic_pipeline/c_trim \
    results/1_basic_pipeline/f_align_and_count
