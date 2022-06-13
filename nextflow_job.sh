#!/bin/bash
#
# This script sets up the nextflow pipeline as a job on
# the imperial computing cluster. Nextflow will run as
# a PBS job, allowing it to download the pipeline from
# github and orchestrate the running of the pipeline
# on PBS. Based on the imperial_rcs config, the pipeline
# results will be saved in a directory named 
# 3_nextflow_pipeline_results.
#
# job specification
#PBS -lselect=1:ncpus=4:mem=16gb
#PBS -lwalltime=06:00:00

# cd to the directory the job was launched from
cd "$PBS_O_WORKDIR"

# activate conda module on the imperial cluster
module load anaconda3/personal

# create conda environment if needed
# conda env remove -n recode_rnaseq
# conda update conda
# conda install mamba
# mamba env create -f environment.yml

# activate conda environment
source activate recode_rnaseq

# export environment variables
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_TEMP="/rds/general/ephemeral/user/$USER/ephemeral/tmp"
export NXF_WORK="/rds/general/ephemeral/user/$USER/ephemeral/tmp/work"
export NXF_SINGULARITY_CACHEDIR="/rds/general/ephemeral/user/$USER/ephemeral/singularity_cache"

# run ImperialCollegeLondon/ReCoDE_rnaseq_pipeline
nextflow run main.nf -profile imperial_rcs
