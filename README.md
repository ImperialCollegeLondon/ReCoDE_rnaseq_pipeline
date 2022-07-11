
ReCoDE RNA-seq
================

<!-- buttons -->
[![actions status](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/workflows/Pipeline%20CI/badge.svg)](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

## Introduction

This repository is part of the REsearch COmputing and Data science Exemplars (ReCoDE) project. ReCoDE exemplars are annotated high-quality research software projects that aim to support teaching and learning in research computing and data science. Courses on these topics usually focus on a single topic at a time, so ReCoDE exemplars help students understand how this knowledge can be integrated into a functional project.

This exemplar involves the development of a pipeline for processing large volumes of biological data. The development of next generation sequencing technologies has facilitated a systems-level approach to biological and biomedical research. In particular, RNA sequencing (RNA-seq) has become a ubiquitous method for gene expression profiling. However, the colossal datasets produced by this method pose a new challenge for life science researchers, who commonly have little statistical or computational training. The processing of sequencing data commonly requires the development of custom workflows that can run in parallel on university computing clusters. Furthermore, the quality control and statistical analysis of these data requires specialist knowledge and purpose-built software packages.

This project demonstrates the development of a pipeline for processing RNA-seq datasets on the computing clusters, such as the Imperial College Research Computing Service (RCS), and basic statistical analysis of the normalised data. This will involve: 

- Quality control and trimming of raw RNA-seq reads
- Alignment of reads to the human reference genome
- Conversion of aligned reads to a matrix of gene counts
- Downstream statistical analysis:
  - Data normalisation 
  - Unsupervised analysis (e.g. PCA)
  - Differential expression and enrichment using edgeR

We will write scripts, using bash and the R programming language, that can execute these steps in parallel on computing clusters. 

![A flow diagram outlining the RNA-seq analysis workflow](assets/flow.png?raw=true "An overview of RNA sequencing, data preprocessing and downstream analysis.")

## Tools used

While there are many tutorials that discuss the processing of RNA-seq data, this exemplar will focus on: i) the parallelisation of bioinformatics tools on computing clusters, such as the Imperial RCS; ii) the development of a reproducible pipeline, introducing methods and tools that can help with this, including:
- [Nextflow](https://www.nextflow.io/), a workflow management system that makes it easy to develop data-driven pipelines.
- [Conda](https://docs.conda.io/en/latest/), an package management system that allows you to share your local environment with others.
- [Docker](https://www.docker.com/), an application for packaging dependencies into a virtual container. 
- [Git/GitHub](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/), a version control system that integrates with nextflow to make pipelines shareable.
- [Continous integration](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions), a practice of automatic code testing.

## Quick start

To get started using the pipeline, install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and either of [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/docs) or [Conda](https://docs.conda.io/en/latest/). Then, the following command can be used to run the pipeline using a small test dataset:

```
nextflow run ImperialCollegeLondon/ReCoDE_rnaseq_pipeline -profile test,docker
```

Note that this pipeline is not designed to handle all types of RNA-seq data (e.g. it was not designed for paired-read data). If you have large amounts of RNA-seq data to process, we recommend using the [nextflow-core RNA-seq pipeline](https://github.com/nf-core/rnaseq).

## Using the repository for learning

This repository was created with in-depth annotations to explain the process of building the pipeline. The documents in [`docs/`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/learning) explain how the data processing pipelines were developed. We created three versions of the pipeline to show how simple pipelines can be integrated and enhanced using Nextflow and other tools. 

1. The first pipeline is a simple bash script that runs each data processing stage sequentially. This method will take a long time to run for larger datasets. See the document [`docs/simple_local_pipeline.md`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/simple_local_pipeline.md) for more information.
2. The second pipeline runs the same steps, however it allows each RNA-seq sample to be run in parallel on the Imperial computing cluster. It is a lot faster, but can be unwieldy to use for large numbers of samples. See the document [`docs/parallelised_pipeline.md`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/parallelised_pipeline.md) for more information. 
3. The full pipeline adapts these steps using Nextflow. This pipeline can be more easily shared with others and will automatically orchestrate the running of large numbers of samples. It will run seamlessly both locally, on computing clusters and on cloud platforms. See the document [`docs/nextflow_pipeline.md`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/nextflow_pipeline.md) for more information.

To download the pipeline, install [Git](https://github.com/git-guides/install-git) and clone this repository from the command line, using the following code:
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

The downstream analysis steps for RNA-seq data require less compute power and often need a more customised workflow. So, this is demonstrated separately in an R markdown notebook. You can run this notebook, located at [`notebooks/downstream_analysis.Rmd`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/notebooks/downstream_analysis.Rmd), or you can view a complete markdown version of the notebook in [`docs/downstream_analysis.md`](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/tree/main/docs/downstream_analysis.md). 

## Reporting issues

If you have any problems running the pipelines, feel free to report them as [an issue](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/issues). If you find any mistakes and know how to correct them, you could also make the modification and create [a pull request](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/pulls).

## License

