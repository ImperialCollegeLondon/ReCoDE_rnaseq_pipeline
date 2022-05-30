#PBS -lselect=1:ncpus=1:mem=4gb
#PBS -lwalltime=00:30:00

# cd to the directory the job was launched from
cd $PBS_O_WORKDIR 

# stop running the script if there are errors
set -e

# where to save the pipeline results
RES_DIR="2_parallelised_pipeline_results"

# load multiqc
module load multiqc/0.6
source activate multiqc

# combine the fastqc results
bin/multiqc.sh \
        $RES_DIR/b_multiqc \
        $RES_DIR/a_fastqc

# combine the fastqc results generated for the trimmed fastq files
bin/multiqc.sh \
        $RES_DIR/d_trim_multiqc \
        $RES_DIR/c_trim

# remove unzipped gtf
rm $RES_DIR/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf

# use multiqc to assess the alignment and counts
bin/multiqc.sh \
    $RES_DIR/g_final_multiqc \
    $RES_DIR/c_trim \
    $RES_DIR/f_align_and_count
