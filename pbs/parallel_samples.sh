#PBS -lselect=1:ncpus=4:mem=32gb
#PBS -lwalltime=01:00:00

# cd to the directory the job was launched from
cd $PBS_O_WORKDIR 

# number of cores available
NUM_CORES=4

# stop running the script if there are errors
set -e

# get names of samples
readarray -t SAMPLE_SRR < data/files.txt

# get name of sample to run
s="${SAMPLE_SRR[$PBS_ARRAY_INDEX - 1]}"

# where to save the pipeline results
RES_DIR="2_parallelised_pipeline_results"

# run fastqc on raw fastq
bin/a_fastqc.sh \
    $RES_DIR/a_fastqc \
    data/fastq/$s.fastq.gz

# trim the fastq files
bin/c_trim.sh \
    $RES_DIR/c_trim \
    data/fastq/$s.fastq.gz \
    $NUM_CORES

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
