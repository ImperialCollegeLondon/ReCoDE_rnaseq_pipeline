#PBS -lselect=1:ncpus=2:mem=32gb
#PBS -lwalltime=02:00:00

# cd to the directory the job was launched from
cd $PBS_O_WORKDIR 

# number of cores available
NUM_CORES=2

# stop running the script if there are errors
set -e

# make sure scripts are executable
chmod u+x data/get_data.sh
chmod u+x bin/*.sh

# get names of samples to run
readarray -t SAMPLE_SRR < data/files.txt

# where to save the pipeline results
RES_DIR="2_parallelised_pipeline_results"

# make the top level results folder
if [ -e  $RES_DIR ]; then
    echo "Results folder already exists, previous files may be overwritten."
else
    mkdir $RES_DIR
fi

# function creates subfolders within the results
create_folder () {
    if [ ! -e  $RES_DIR/$1 ]; then
        mkdir $RES_DIR/$1
    fi
}

create_folder "a_fastqc"
create_folder "c_trim"
create_folder "e_star_index"

# load STAR on the cluster
module load star/2.7.1a

# index the genome using STAR
bin/star_index.sh \
        $RES_DIR/e_star_index/ \
        data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz \
        data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz \
        1 \
        $NUM_CORES

# remove unzipped fasta
rm $RES_DIR/e_star_index/GCF_000004515.6_Glycine_max_v4.0_genomic.fna

create_folder "f_align_and_count"
