#!/bin/bash
#
# This script is run for each sample in parallel on the cluster by the PBS
# QC is performed, reads are trimmed before being mapped and quantified
#
# job specification
#PBS -lselect=1:ncpus=2:mem=16gb
#PBS -lwalltime=01:00:00

# set to data for the full dataset or $DATA_DIR for the test dataset
DATA_DIR="data/test"

# cd to the directory the job was launched from
cd "$PBS_O_WORKDIR"

# activate conda module on the imperial cluster
module load anaconda3/personal

# activate conda environment
source activate recode_rnaseq

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

# run fastqc on raw fastq
bin/fastqc.sh \
  "${RES_DIR}/a_fastqc" \
  "$DATA_DIR/fastq/${s}.fastq.gz" \
  "${RES_DIR}/a_fastqc/${s}"

# trim the fastq files
bin/trim.sh \
  "${RES_DIR}/c_trim" \
  "${DATA_DIR}/fastq/${s}.fastq.gz" \
  "${RES_DIR}/c_trim/${s}"

# perform alignment using STAR, providing the directory of the indexed genome
bin/align.sh \
  "${RES_DIR}/e_star_index" \
  "${RES_DIR}/c_trim/${s}_trimmed.fq.gz" \
  "${RES_DIR}/f_align/${s}" \
  "${NUM_CORES}"

# remove unzipped fastq
rm "${RES_DIR}/f_align/${s}.fastq"

bin/count.sh \
  "${RES_DIR}/f_align/${s}Aligned.sortedByCoord.out.bam" \
  "${RES_DIR}/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf" \
  "${RES_DIR}/g_count/${s}"

# check the counts have been successfully created
if [ -e  "${RES_DIR}/g_count/${s}.counts" ]; then
  echo "Counts file for sample ${s} was successfully created"
else
  err "Counts file for sample ${s} was not created"
fi
