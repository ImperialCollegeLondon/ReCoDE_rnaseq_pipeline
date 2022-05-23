
# a vector of the sequence read archive identifiers for the RNA-seq samples
#            1           2           3           4           5           6
SAMPLE_SRR=("SRR391535" "SRR391536" "SRR391537" "SRR391538" "SRR391539" "SRR391541")

# load fastq-dump command on the cluster
# module load sra-toolkit/2.8.1

# create a file to list the IDs and file locations
rm -f data/files.txt
touch data/files.txt

# for each sample
for srr in ${SAMPLE_SRR[@]}; do

    # print the sample name
    echo $srr

    # create file listing the ID and file location
    echo $srr >> data/files.txt

    # get the sample's file from the sequence read archive
    fastq-dump --gzip -O data/fastq/ "$srr"

done

# get the genome and annotation files from NCBI
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000004515.6

mkdir data/genome

wget -O data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz \
     "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz"

wget -O data/genome/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz \
     "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.gtf.gz"
