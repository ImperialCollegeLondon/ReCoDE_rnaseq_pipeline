
# if necessary, load the relevant command line tool on the computing cluster
# module load sratoolkit/2.8.1
# export PATH=$PATH:~/apps/sratoolkit.3.0.0-ubuntu64/bin

# a vector of the SRR identifiers for the RNA-seq samples
#            1           2           3           4           5           6
SAMPLE_SRR=("SRR391535" "SRR391536" "SRR391537" "SRR391538" "SRR391539" "SRR391541")

# create a file to list the IDs and file locations
touch data/files.txt

# for each sample
for srr in ${SAMPLE_SRR[@]}; do

    # print the sample name
    echo $srr

    # create file listing the ID and file location
    echo $srr >> data/files.txt

    # get the sample's file
    fastq-dump --gzip "$srr"

done
