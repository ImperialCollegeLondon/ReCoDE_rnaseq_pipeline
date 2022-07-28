#!/bin/bash
#
# This script sets up the jobs in pbs/ on the imperial cluster
# they are set as dependent jobs, so they should run one after the other
# if you want to change the data to be run, change the data_dir variable in pbs/parallel_samples.pbs and setup_pipeline.pbs

# setup the pipeline - including getting the data, setting up folders and indexing the genome
sp_jid="$(qsub pbs/setup_pipeline.pbs)"

# now we have the data, get a list of the samples
while IFS=\= read srr; do
    SAMPLE_SRR+=($srr)
done < data/files.txt

# run the pipeline in parallel for each sample
ps_jid="$(qsub -W depend=afterok:"${sp_jid}" -J 1-"${#SAMPLE_SRR[@]}" pbs/parallel_samples.pbs)"

# use multiqc to group the QC files
mq_jid="$(qsub -W depend=afterok:"${ps_jid}" pbs/multiqc.pbs)"
