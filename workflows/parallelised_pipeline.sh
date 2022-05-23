
# setup the pipeline - including getting the data, setting up folders and indexing the genome
sp_jid=$(qsub pbs/setup_pipeline.sh)

# now we have the data, get a list of the samples
readarray -t SAMPLE_SRR < data/files.txt

# run the pipeline in parallel for each sample
ps_jid=$(qsub -W depend=afterok:$sp_jid -J 1-${#SAMPLE_SRR[@]} pbs/parallel_samples.sh)

# use multiqc to group the QC files
mq_jid=$(qsub -W depend=afterok:$ps_jid pbs/multiqc.sh)
