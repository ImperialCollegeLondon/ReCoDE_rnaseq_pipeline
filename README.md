
RNA-seq Analysis
================

<!-- buttons -->
[![actions status](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/workflows/Pipeline%20CI/badge.svg)](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

### Author: Jack Gisby

## Description

The RNA-seq analysis exemplar involves the development of a pipeline for processing large volumes of biological data. The development of next generation sequencing technologies has facilitated a systems-level approach to biological and biomedical research. In particular, RNA sequencing (RNA-seq) has become a ubiquitous method for gene expression profiling. However, the colossal datasets produced by this method pose a new challenge for life science researchers, who commonly have little statistical or computational training. The processing of sequencing data commonly requires the development of custom workflows that can run in parallel on university computing clusters. Furthermore, the quality control and statistical analysis of these data requires specialist knowledge and purpose-built software packages.

This project demonstrates the development of a pipeline for processing RNA-seq datasets on the computing clusters, such as the Imperial College Research Computing Service (RCS), and basic statistical analysis of the normalised data. This will involve: 

- Quality control and trimming of raw RNA-seq reads
- Alignment of reads to the human reference genome
- Conversion of aligned reads to a matrix of gene counts
- Downstream statistical analysis:
  - Data normalisation 
  - Unsupervised analysis (e.g. PCA)
  - Differential expression and enrichment using edgeR


## Prerequisites:

- Familiarity with bash (A course such as "[The Linux Command Line for Scientific Computing](https://www.imperial.ac.uk/study/pg/graduate-school/students/doctoral/professional-development/research-computing-data-science/courses/linux-command-line-for-scientific-computing/)", hosted by the Imperial Research Computing & Data Science Team, would be provide a suitable background.)
- Familiarity with R programming language
- Familiarity with computing clusters 

## Learning Outcomes:

Upon completion of this tutorial, students will be able to:

1. parallelise bioinformatics tools on computing clusters, such as the Imperial RCS
2. develop a reproducible pipeline

Tools used achieve this as part of the exemplar include:
- [Nextflow](https://www.nextflow.io/), a workflow management system that makes it easy to develop data-driven pipelines.
- [Conda](https://docs.conda.io/en/latest/), an package management system that allows you to share your local environment with others.
- [Docker](https://www.docker.com/), an application for packaging dependencies into a virtual container. 
- [Git/GitHub](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/), a version control system that integrates with nextflow to make pipelines shareable.
- [Continous integration](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions), a practice of automatic code testing.


## Project Structure

In this exemplar, we set up three versions of the same pipeline, that process RNA sequencing data. They are available in [docs](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/), along with the background information necessary to follow each pipeline. Each of these pipeline generates .count files in their final step, which counts how many RNA sequences map to each gene for each sample. We created multiple versions of the pipeline to show the different ways in which researchers might want to process their data. They are as follows:

1. [simple_pipeline](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/simple_local_pipeline.md): This is the simplest version of the pipeline. It is a bash script that can be run straight from the command line. It will execute each stage of the pipeline in sequentially. 

    Pros:

    * Simple to follow and execute

    Cons:

    * slow for a large number of samples
    * if there is an error for a single sample the entire pipeline will fail


2. [parallelised_pipeline](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/parallelised_pipeline.md): This builds on the simple pipeline to run the pipeline on the cluster. It still uses bash scripts to orchestrate the pipeline, however it also allows each sample to be run in parallel on the cluster. This version of the pipeline is therefore a lot faster!

    Pros:

    * a lot faster than simple_pipeline

    Cons:

    * unwieldy to use for large numbers of samples    


3. [nextflow_pipeline](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/nextflow_pipeline.md): This is the most advanced version of the pipeline. Instead of relying on bash scripts, that could fail without warning, it uses the nextflow workflow manager to run the pipeline stages either on your local computer or on a computing cluster. 

    Pros:

    * automatically orchestrate the running of large numbers of samples
    * run seamlessly both locally, on computing clusters and on cloud platforms
    * easily shared with others

    Cons:

    * takes more time to set up than the simple or parallelised pipeline 

We created markdown documents, available in this directory, for each of these pipelines. These documents explain the pipelines and link to external resources in case you want to learn more. We suggest you go through the documentation in the order above, as each one builds upon the last.

The downstream analysis uses the output of these pipelines to investigate the data. The dataset is a lot smaller at this point, so this stage of the analysis can be accomplished locally without the use of the cluster. The code for this stage is contained within the notebooks/ directory. The code is written in an Rmarkdown document, which has been run and stored as a markdown document within this directory. 



## Getting started

To get started using the pipeline, install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and either of [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/docs) or [Conda](https://docs.conda.io/en/latest/). Then, the following command can be used to run the pipeline locally using a small test dataset:

```
nextflow run -r main ImperialCollegeLondon/ReCoDE_rnaseq_pipeline -profile test,docker
```

Note that this pipeline is not designed to handle all types of RNA-seq data (e.g. it was not designed for paired-read data). If you have large amounts of RNA-seq data to process, we recommend using the [nextflow-core RNA-seq pipeline](https://github.com/nf-core/rnaseq).

To download the pipelines, install [Git](https://github.com/git-guides/install-git) and clone this repository from the command line, using the following code:
```
git clone https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline.git
```

To run the first two pipelines, you will need to install the command line applications listed in [`environment.yml`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/blob/main/environment.yml). Alternatively, you can use conda to install these applications, using the following code:
```
conda env create -f environment.yml
conda activate recode_rnaseq
```

This conda environment also contains Nextflow, which we will be using to create the final iteration of the data processing pipeline. To run the nextflow pipeline, you will also need either Docker or Conda installed. 

The scripts in this repository are designed for use with the Imperial computing cluster. If you do not have access to the cluster, you might be able to adapt the code to your own cluster's configuration. Alternatively, the primary pipeline uses nextflow, which is adaptable to many different platforms. You could run the nextflow pipeline on your local computer, or configure it to run on another cluster or even the cloud.

Some of these tools are not available on windows. If you wish to run the basic pipeline discussed in this document but you use Windows, you could use the Windows Subsystem for Linux for the course or work on a computing cluster, such as [Imperial's high performance computing cluster](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/support/getting-started/
). 

The downstream analysis steps for RNA-seq data require less compute power and often need a more customised workflow. So, this is demonstrated separately in an R markdown notebook. You can run this notebook, located at [`notebooks/downstream_analysis.Rmd`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/notebooks/downstream_analysis.Rmd), or you can view a complete markdown version of the notebook in [`docs/downstream_analysis.md`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/downstream_analysis.md). 

## Reporting issues

If you have any problems running the pipelines, feel free to report them as [an issue](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/issues). If you find any mistakes and know how to correct them, you could also make the modification and create [a pull request](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/pulls).

## License

