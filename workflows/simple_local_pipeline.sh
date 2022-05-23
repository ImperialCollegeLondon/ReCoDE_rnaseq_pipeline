
# stop running the script if there are errors
set -e

# download the data (skip this step if already downloaded)
if [ ! -f data/files.txt ]; then
    data/get_data.sh
fi

# get names of samples to run
readarray -t SAMPLE_SRR < data/files.txt

# where to save the pipeline results
RES_DIR="1_simple_local_pipeline_results"

# number of cores available
NUM_CORES=6

# make the top level results folder
if [ -e  $RES_DIR ]; then
    echo "Results folder already exists, previous files may be overwritten."
else
    mkdir $RES_DIR
fi

# function creates subfolders within the results
create_folder () {
    if [ ! -e  $RES_DIR/$1 ]; then
        mkdir $RES_DIR/$1
    fi
}

# make a subfolder for the fastqc output
create_folder "a_fastqc"
create_folder "c_trim"

# for each sample
for s in "${SAMPLE_SRR[@]}"; do

    # run fastqc on raw fastq
    bin/a_fastqc.sh \
        $RES_DIR/a_fastqc \
        data/fastq/$s.fastq.gz \
        $s
done

# combine the fastqc results
bin/b_multiqc.sh \
        $RES_DIR/b_multiqc \
        $RES_DIR/a_fastqc

for s in "${SAMPLE_SRR[@]}"; do

    # trim the fastq files
    bin/c_trim.sh \
        $RES_DIR/c_trim \
        data/fastq/$s.fastq.gz \
        $s
done

# combine the fastqc results generated for the trimmed fastq files
bin/d_trim_multiqc.sh \
        $RES_DIR/d_trim_multiqc \
        $RES_DIR/c_trim

create_folder "e_star_index"

# index the genome using STAR
bin/e_star_index.sh \
        $RES_DIR/e_star_index/ \
        data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz \
        data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz \
        3 \
        $NUM_CORES

# remove unzipped fasta
rm $RES_DIR/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.fna

create_folder "f_align_and_count"

for s in "${SAMPLE_SRR[@]}"; do

    # perform alignment using STAR, providing the directory of the indexed genome
    bin/f_align_and_count.sh \
        $RES_DIR/e_star_index \
        $RES_DIR/c_trim/"$s"_trimmed.fq.gz \
        $RES_DIR/f_align_and_count/$s \
        $NUM_CORES \
        $RES_DIR/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

    # remove unzipped fastq
    rm $RES_DIR/f_align_and_count/$s.fastq 

    # check the counts have been successfully created
    if [ -e  $RES_DIR/f_align_and_count/$s.counts ]; then
        success=""
    else
        success=" not"
    fi

    echo "Counts file for sample $s was$success successfully created"

    # TODO: check counts file is non-empty
done

# remove unzipped gtf
rm $RES_DIR/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

# use multiqc to assess the alignment and counts
bin/g_final_multiqc.sh \
    $RES_DIR/g_final_multiqc \
    $RES_DIR/c_trim \
    $RES_DIR/f_align_and_count
