
readarray -t SAMPLE_SRR < data/files.txt

sp_jid=$(qsub pbs/setup_pipeline.sh)
ps_jid=$(qsub -W depend=afterok:$sp_jid -J ${#SAMPLE_SRR[@]} parallel_samples.sh)
mq_jid=$(qsub -W depend=afterok:$ps_jid multiqc.sh)
