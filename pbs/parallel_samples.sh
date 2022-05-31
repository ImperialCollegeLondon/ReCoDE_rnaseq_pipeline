#PBS -lselect=1:ncpus=2:mem=16gb
#PBS -lwalltime=01:00:00

# cd to the directory the job was launched from
cd $PBS_O_WORKDIR 

# number of cores available
NUM_CORES=2

# stop running the script if there are errors
set -e

# get names of samples
readarray -t SAMPLE_SRR < data/files.txt

# get name of sample to run
s="${SAMPLE_SRR[$PBS_ARRAY_INDEX - 1]}"

# where to save the pipeline results
RES_DIR="2_parallelised_pipeline_results"

# load fastqc
module load fastqc/0.11.9

# run fastqc on raw fastq
bin/fastqc.sh \
    $RES_DIR/a_fastqc \
    data/fastq/$s.fastq.gz \
    $RES_DIR/a_fastqc/$s

# load trimgalore
module load trim_galore/0.4.1
module load cutadapt/1.9.1

# trim the fastq files
bin/trim.sh \
    $RES_DIR/c_trim \
    data/fastq/$s.fastq.gz \
    $RES_DIR/c_trim/$s

# load STAR and htseq
module load star/2.7.1a 
module load htseq/0.6.1

# perform alignment using STAR, providing the directory of the indexed genome
bin/align_and_count.sh \
    $RES_DIR/e_star_index \
    $RES_DIR/c_trim/"$s"_trimmed.fq.gz \
    $RES_DIR/f_align/$s \
    $NUM_CORES

# remove unzipped fastq
rm $RES_DIR/f_align/$s.fastq 

bin/count.sh \
    $RES_DIR/f_align/"$s"Aligned.sortedByCoord.out.bam \
    $RES_DIR/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf \
    $RES_DIR/g_count/$s 

# check the counts have been successfully created
if [ -e  $RES_DIR/f_align/$s.counts ]; then
    success=""
else
    success=" not"
fi

echo "Counts file for sample $s was$success successfully created"
