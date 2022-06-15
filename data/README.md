# Soybean data

The data we will be using were deposited in the Sequence Read Archive (SRA) under the identifier SRP009826, entitled "Comparing the response to chronic ozone of five agriculturally important legume species."

RNA-seq was performed for soybean plants that were part of the control group (ambient O3 levels) and treatment group (elevated O3 levels).

This data was used previously in "The bench scientist's guide to statistical analysis
of RNA-Seq data" (http://www.biomedcentral.com/1756-0500/5/506) and a tutorial (https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/).

RNA-seq data can be quite large. Each sample might generate a few gigabytes of data. So, if you perform an experiment with hundreds of samples, your dataset might be in the order of terabytes! If you wish to try and process the full RNA-seq data, you can use the script `data/get_data.sh` to download the full data. However, these data may be too large to run on your laptop! For this reason, we have taken a subset of the full soybean data, which is included in this repository. These data are available in the `data/test/` directory.

The file 'get_data.sh' is used to download the data from the SRA. We use the command line tool `fastq-dump` and provide the identifiers of the individual samples we want to download.

