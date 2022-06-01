
#!/bin/bash
#
# This script sets up the jobs in pbs/ on the imperial cluster
# they are set as dependent jobs, so they should run one after the other
# if you want to change the data to be run, change the data_dir variable in pbs/parallel_samples.sh and setup_pipeline.sh

# setup the pipeline - including getting the data, setting up folders and indexing the genome
sp_jid="$(qsub pbs/setup_pipeline.sh)"

# now we have the data, get a list of the samples
readarray -t SAMPLE_SRR < data/files.txt

# run the pipeline in parallel for each sample
ps_jid="$(qsub -W depend=afterok:"${sp_jid}" -J 1-"${#SAMPLE_SRR[@]}" pbs/parallel_samples.sh)"

# use multiqc to group the QC files
mq_jid="$(qsub -W depend=afterok:"${ps_jid}" pbs/multiqc.sh)"
