# Learning

This directory contains markdown documents explaining the pipeline and downstream analysis. In this exemplar, we set up three versions of the same pipeline, that process RNA sequencing data. Each of these pipeline generates .count files in their final step, which counts how many RNA sequences map to each gene for each sample. We created multiple versions of the pipeline to show the different ways in which researchers might want to process their data. They are as follows:

 - **simple_pipeline**: This is the simplest version of the pipeline. It is a bash script that can be run straight from the command line. It will execute each stage of the pipeline in sequentially. While this pipeline is very simple, it will be slow for a large number of samples. Also, if there is an error for a single sample the entire pipeline will fail.

 - **parallelised_pipeline**: This builds on the simple pipeline to run the pipeline on the cluster. It still uses bash scripts to orchestrate the pipeline, however it also allows each sample to be run in parallel on the cluster. This version of the pipeline is therefore a lot faster!

 - **nextflow_pipeline**: This is the most advanced version of the pipeline. Instead of relying on bash scripts, that could fail without warning, it uses the nextflow workflow manager to run the pipeline stages either on your local computer or on a computing cluster. 

 We created markdown documents, available in this directory, for each of these pipelines. These documents explain the pipelines and link to external resources in case you want to learn more. We suggest you go through the documentation in the order above, as each one builds upon the last.

The downstream analysis uses the output of these pipelines to investigate the data. The dataset is a lot smaller at this point, so this stage of the analysis can be accomplished locally without the use of the cluster. The code for this stage is contained within the notebooks/ directory. The code is written in an Rmarkdown document, which has been run and stored as a markdown document within this directory. 
