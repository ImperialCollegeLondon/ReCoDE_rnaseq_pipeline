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

# activate conda module on the imperial cluster
module load anaconda3/personal

# create conda environment if needed
# conda env remove -n recode_rnaseq
# conda update conda
# conda install mamba
# mamba env create -f environment.yml

# activate conda environment
source activate recode_rnaseq

nextflow ImperialCollegeLondon/ReCoDE_rnaseq_pipeline -profile test,imperial_rcs
