<!-- buttons -->
[![actions status](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/workflows/test_pipeline/badge.svg)](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

.. |Build Status| image:: https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/workflows/test_pipeline/badge.svg
   :target: https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions
# A cluster-based pipeline for the analysis of RNA-sequencing data

The development of next generation sequencing technologies has facilitated a systems-level approach to biological and biomedical research. In particular, RNA sequencing (RNA-seq) has become a ubiquitous method for gene expression profiling. However, the colossal datasets produced by this method pose a new challenge for life science researchers, who commonly have little statistical or computational training. The processing of sequencing data commonly requires the development of custom workflows that can run in parallel on university computing clusters. Furthermore, the quality control and statistical analysis of these data requires specialist knowledge and purpose-built software packages.

In this exemplar, we will demonstrate the development of a pipeline for processing RNA-seq datasets on the Imperial College Research Computing Service (RCS) and basic statistical analysis of the normalised data. This will involve: 

- Quality control and trimming of raw RNA-seq reads
- Alignment of reads to the human reference genome
- Conversion of aligned reads to a matrix of gene counts
- Downstream statistical analysis:
  - Data normalisation 
  - Unsupervised analysis (e.g. PCA)
  - Differential expression and enrichment using edgeR

We will write scripts, using bash and the R programming language, that can execute these steps in parallel on the RCS. 

![A flow diagram outlining the RNA-seq analysis workflow](assets/flow.png?raw=true "An overview of RNA sequencing, data preprocessing and downstream analysis.")

While there are many tutorials that discuss the processing of RNA-seq data, this exemplar will focus on: i) the parallelisation of bioinformatics tools on the Imperial RCS; ii) the development of a reproducible pipeline, demonstrating the use of environment management systems like conda and potentially introducing workflow languages like nextflow.  

