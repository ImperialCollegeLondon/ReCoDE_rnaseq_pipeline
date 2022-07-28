# Running the parallelised pipeline on Imperial computing cluster

In `docs/simple_local_pipeline.md`, we discussed the limitations of running our pipeline as a simple script on a local machine. We noted that running the pipeline on a standard laptop was only feasible for small datasets, because larger datasets would require more memory than we had available. We could fix this by running the pipeline on the cluster instead as a job - this way, we can access as much memory as we may need. So, why do we need to modify the pipeline, beyond just adding resource specifications to the top of the script?

We also said that the pipeline was limited because it could only run a single task at a time, so it will be slow when we have large numbers of samples. This is a more difficult problem to solve. Some of our pipeline steps, for instance indexing the genome with STAR and generating reports with MultiQC, only need to be run once. Other steps, like FastQC and STAR alignment, must be run for every single sample. So, we will have to reconfigure our pipeline to run efficiently on the cluster.

## Set up

To run the full pipeline on the cluster, you must have setup the conda environment and downloaded this repository. The instructions to do this are available in this exemplar's top-level `README.md`. Note, however, that conda works slightly differently on the imperial cluster. You must activate the conda module using the following code:
```
module load anaconda3/personal
```

If you have never run conda before, you must initialise it like so:
```
anaconda-setup
```

Finally, you can install the recode_rnaseq environment:
```
conda update conda
conda env create -f environment.yml
```

Note that to run the pipeline and to install the conda environment, you must run the commands from within the ReCoDE_rnaseq_pipeline repository. The code to setup the conda environment is also included in the script `pbs/setup_pipeline.pbs`, so you could just uncomment those lines and run the pipeline and install the environment in one go. The install commands in this script use `mamba`, which is a lot faster than conda. So, if you're having trouble installing the environment, try using `mamba`. 

You can then run the pipeline either by running the command `workflows/parallelised_pipeline.sh` from the cluster command line (do not use the qsub command to run this script), or by running the individual lines of code within `workflows/parallelised_pipeline.sh`. On the Imperial cluster, this should automatically run the pipeline for the test dataset. If you are using a cluster owned by a different institution, you may need to do some work to adapt the scripts to your scheduler. If you have any problems, please raise them as [an issue](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/issues).

## Adapting the simple pipeline to run on the cluster

Our first step is to break the pipeline up into multiple stages. We need to separate the elements of the pipeline that should only be run once, and the parts of the pipeline that need to be run for each sample. Specifically, we need to create the results folders before we start generating results. We also need to create a STAR index before we run the alignment step. So, we create a script, `pbs/setup_pipeline.pbs`, which runs these stages first. While MultiQC is another step that only needs to be run once, it can't be run in this script because it has to be run last! 

So, the next step is to carry out the data processing steps for each sample. In `workflows/simple_local_pipeline.sh`, we used loops to run the QC, trimming, alignment and counting scripts. Since we are now on the cluster, we will instead leverage array jobs to run these steps for each of our samples. We create a script, `pbs/parallel_samples.pbs`, which is designed to process a single sample. 

In the `workflows/simple_local_pipeline.sh` script, we looped through the sample IDs, and in each iteration we assigned the sample ID to the variable `s`. In `pbs/parallel_samples.pbs`, we take all the code from within these loops, but we instead define `s` at the top of the script using the environment variable `PBS_ARRAY_INDEX`. We will later set up this script to run as an array, such that `pbs/parallel_samples.pbs` is run six times, once for each sample; each of the six array jobs will have a different number for the variable `PBS_ARRAY_INDEX`, so that each script has a different value of `s`. To do this, we get the list of samples names like we did in the simple pipeline, as follows:
```
while IFS=\= read srr; do
    SAMPLE_SRR+=($srr)
done < data/files.txt
```

Then, we select the Nth value of the sample list (i.e. the Nth row of `data/files.txt`), in the following code:
```
s="${SAMPLE_SRR[$PBS_ARRAY_INDEX - 1]}"
```

The script then runs the data processing steps for the corresponding sample. Once all of the samples have been processed, we can finally run MultiQC to collect the reports. This is achieved using the script `pbs.multiqc.pbs`. 

We also added the PBS parameters to the top of each of the `.pbs` scripts. We gave both `pbs/setup_pipeline.pbs` and `parallel_samples.pbs` multiple cores and high amounts of memory, because STAR can make use of them. MultiQC is a lot less demanding, so we gave `pbs/multiqc.pbs` only a single core and 4GB of RAM. These parameters might need modifications for different genomes or larger amounts of data, but for the soybean dataset they are sufficient. 

## Running the pipeline

Having now split `workflows/simple_local_pipeline.sh` into three different stages, we can submit them to the scheduler. The code needed to do this is available in `workflows/parallelised_pipeline.sh`. This contains the `qsub` commands we need to submit each of the job to the cluster in order. The first command in this script, below, submits the setup script to the scheduler and stores its job ID to the variable `sp_jid`. 
```
sp_jid="$(qsub pbs/setup_pipeline.pbs)"
```

We next read the samples to be processed as a bash array, as below. We do this to calculate the number of samples to be processed, so that we know how large the array job needs to be. 
```
while IFS=\= read srr; do
    SAMPLE_SRR+=($srr)
done < data/files.txt
```

The most complex command is the submission of the array job:
```
ps_jid="$(qsub -W depend=afterok:"${sp_jid}" -J 1-"${#SAMPLE_SRR[@]}" pbs/parallel_samples.pbs)"
```

Like we did for the setup script, we submit `pbs/parallel_samples.pbs` to the scheduler using `qsub` and save its job ID to the variable `ps_jid`. Conveniently, all of the array jobs are connected to a single job ID, which we can store to track its progress. The argument `-W depend=afterok:"${sp_jid}"` tells the scheduler not to start the array job until the setup job is complete. Finally, the `1-"${#SAMPLE_SRR[@]}"` argument sets the script up as an array job with array IDs from 1 to the number of samples in the `SAMPLE_SRR` array. 

Finally, we schedule the MultiQC script to run after the array job is complete, like so:
```
mq_jid="$(qsub -W depend=afterok:"${ps_jid}" pbs/multiqc.pbs)"
```

## Running the pipeline for the full dataset

Now we can try and run the pipeline for the full soybean dataset! This is a lot bigger, but now that we are on the cluster it should be simple. The first step is to change the `DATA_DIR` variable in each of the scripts in the `pbs/` directory from "data/test" to "data/". Then, you will need to download the data! 

There is a script available for doing this: `data/get_data.sh`. This will take some time, so you will need to run it as a PBS job. The easiest way to do this is to uncomment the code `data/get_data.sh` in the `pbs/setup_pipeline.pbs`. This way you will download the data as part of the setup - you may need to increase the `lwalltime` parameter so that everything runs though! Note that this script uses a tool called `fastq-dump` to get the data from the Sequence Read Archive. You could either try installing it using conda, or on the imperial cluster you can use the sra-toolkit module. To do this, just uncomment the line `module load sra-toolkit/2.8.1` in the `data/get_data.sh` script. 

Now that you have done this, you can run the script like we did for the test data. It will take a bit longer, but not too long because we are only running six samples. The benefit of the parallelised pipeline would be even more obvious if we were to try running hundreds of samples! Once the pipeline has run, you can check the MultiQC outputs and see how they look for a real dataset. If of interest, you could try creating a job to run `workflows/simple_local_pipeline.sh` for the full dataset. Since we only have six samples, this will run, but it will be a lot slower than the parallelised pipeline!

## Limitations

The parallelised pipeline alleviated some of the limitations of the simple local pipeline. Running the pipeline on the cluster gives us access to the resources we need to run compute-intensive tools like STAR. We also made use of array jobs to run each sample in parallel, making the processing a lot faster. 

In `docs/simple_local_pipeline.md` we also noted that the pipeline was prone to failure. What if a file wasn't generated for one of the samples and we never noticed? What if a file was generated but there was an error that prevented the process from finishing? These problems may have been missed by the simple pipeline. And, sadly, they may also be missed by the parallelised pipeline. If the array job fails completely, the job will finish and MultiQC won't run, because it depends on the array job. So, we might pick up major problems with the pipeline. But, there might still be an error for a particular sample that doesn't get picked up. In this case, the parallelised pipeline might be even worse than the simpler one! To make sure each stage ran correctly we would have to check each of the error files produced by each of the array jobs. 

We could make this pipeline more robust. For instance, we could improve the final script (`pbs/multiqc.pbs`) by adding code that checks the right files were generated for all of the samples. This would make us a lot more confident that the pipeline had run correctly. But, writing all these checks takes time that we could be using to process more data or perform downstream analysis. And, eventually we will write a pipeline and forget to check something important. 

In the next document, `docs/nextflow_pipeline.sh`, we will discuss the use of Nextflow to manage our pipeline. Nextflow will allow us to more elegantly tie our pipeline steps together, and it will make sure that every stage produces the expected outputs! 
