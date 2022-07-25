# Computing clusters

Often, as is often the case with our RNA-seq data, we need more computing power than is available on a typical laptop. RNA-seq data can take an extremely long time to process on a single machine. Furthermore, lots of memory (RAM) is required for running alignment tools. You might also find that you just don't have enough space to store a large dataset on your computer! 

For this reason, many universities have computing clusters that can help overcome these issues. Computing clusters are a set of computers that are linked together by a network so that they can work together. Each of the computers in the cluster is referred to as a node. These clusters are usually set up with some sort of job scheduler, which allows researchers to submit code to run as a "job"; the scheduler will find a slot on one or more of the nodes on which the code can be run. The Imperial computing cluster uses the Portable Batch System (PBS) as a scheduler. 

In this document, we will demonstrate how the RNA-seq pipeline can be set up to run on the Imperial computing cluster. We will expect that you will have had some previous exposure to running jobs on the cluster, for instance by taking the course [Introduction to HPC at Imperial](https://www.imperial.ac.uk/study/pg/graduate-school/students/doctoral/professional-development/research-computing-data-science/courses/introduction-to-hpc/) hosted by the Research Computing and Data Science team. If you have access to a cluster other than the Imperial cluster, you could try adapting the scripts for compatibility with your cluster's scheduler. If you do not have access to a computing cluster, it might be difficult to follow along with this document; we suggest trying to understand the strengths and limitations of setting the pipeline up on the cluster, then moving on to `docs/nextflow_pipeline.md` which describes the final iteration of our pipeline, and can be run without access to a computing cluster.

## How to run scripts on the cluster

The Imperial cluster uses the PBS scheduler. A job can be started by running the `qsub` command to queue a script to run, for instance as follows:
```
qsub my_script.pbs
```

The script is setup in a `.pbs` file. This script is very similar to a standard bash script, however it must contain some extra information that tells the scheduler what computing resources the script will need. The key settings you will need to change are the length of time the job will require, the number of cores that should be allocated to the job and the total memory the job will need. For instance, if you were to run a job for 2 hours, with 3 cores and 4GB of memory, you could do so by including the following lines near the top of your `.pbs` file:
```
#PBS -lselect=1:ncpus=3:mem=4gb
#PBS -lwalltime=02:00:00
```

The other parameter (`lselect`) is set to 1, indicating the job will run on a single node. It is usually best to stick with a single node, unless you are using an application that can make use of multiple. 

Note that there are many more options we can send to the scheduler, these are just a few that we will commonly want to manipulate. Also, we can pass parameters to PBS in multiple ways. We can include them in our script, as we have done above. This is convenient because we can store our parameters with our script. We know the parameters we used to run the pipeline, and we can use them again next time. Alternatively, these parameters can be passed as arguments to `qsub` using the command line. This can be useful sometimes, but you might want to save a log of the commands you used somewhere!

There are difference "resources" available on the computing cluster. For instance, at time of writing, the "throughput72" resource on the Imperial computing cluster is applicable for jobs that run on a single node, use 1-8 cores, use 1-100GB memory and run for 9-72 hours. Generally, you should pick the smallest runtime, cores and memory that will allow your code to run, because this will lead to the shortest queue time. Most clusters have a variety of resources available. Some may be optimised for smaller jobs; others may be setup to accommodate high memory requirements or the use of many cores; some will allow access to additional resources, such as GPUs. 

The above parameters are great for running a single task on the cluster. Sometimes, however, we might have a sequence of jobs that we want to run, one after the other. Say we want to run three programs, one after the other, named A and B. We could create a `.pbs` script that runs these tasks sequentially, then submit this to the PBS scheduler. However, A and B might have very different resource requirements. Perhaps A requires lots of cores to run in a reasonable amount of time, while task B requires less processing but needs a lot of memory. For this use case, we can set up job dependencies to tell the PBS scheduler to start one job after the other completes. We would first create a `.pbs` script to run programs A (`a.pbs`) and B  (`b.pbs`) with their specific resource requirements, then we could run the following commands:
```
jid_a="$(qsub a.pbs)"

qsub -W depend=afterok:"${jid_a}" b.pbs
```

The first line includes the code `qsub a.pbs`, which schedules `a.pbs` to run like we have seen previously. Running `qsub a.pbs` returns the ID of the job that has been queued. The remaining code on the first line captures this ID into the variable `jid_a`. 

The second line submits `b.pbs` to the scheduler. We add `-W depend=afterok:"${jid_a}"` to ask the scheduler not to start this second job until the first job, with ID `jid_a`, has completed successfully. 

Job dependencies are useful, but could get complex for large numbers of jobs. For instance, in the case of our pipeline, we might want to run the same job for each of our RNA-seq samples. We could have hundreds of samples, meaning we have to manually submit as many jobs! The scheduler's answer to this use case comes in the form of array jobs. Array jobs are created using the J parameter, like so: `-J 1-N`. The scheduler will submit N jobs; these jobs are otherwise identical, except they each have access to the array ID through the environment variable `PBS_ARRAY_INDEX`. This variable will be equal to 1 for job 1 and equal to N for job N. Therefore, we can setup a script that uses the array index to process the Nth RNA-seq sample. Array jobs are convenient because they package a large number of tasks into a single job. We will demonstrate its use for processing RNA-seq data in a later section.

# Environment management

conda lets us create isolated environments that include the packages we need and their dependencies. conda lets us more easily share our code with others, because we can share our code alongside a list of software packages that are required. The easiest way to do this is with an `environment.yml` file, like the one we have created for this exemplar:
```
name: recode_rnaseq
channels:
- bioconda
- conda-forge
dependencies:
- star=2.7.10a
- fastqc=0.11.9
- htseq=0.11.3
- multiqc=1.12
- trim-galore=0.6.7
- sra-tools=2.11.0
- nextflow=22.04.0
- samtools=1.6
```

If we run the code `conda create env -f environment.yml`, it creates an environment with the name "recode_rnaseq". To create this environment, conda installs all the packages listed as dependencies and makes sure that all of the package versions are compatible with each other. conda knows where these packages are located because we supply the channels `bioconda`, a repository of bioinformatics software packages, and `conda-forge`, another repository for open source packages.

conda has some limitations though. In our case, we want to create a reproducible pipeline that can be run easily by anyone. For instance, we might find that less common packages are not available on conda. Or, we might find that the packages we want are only available on certain operating systems; STAR, for instance, cannot run on Windows. We can avoid some of these limitations with docker, that is also well-integrated with nextflow. 

conda allows us to install a set of packages into an environment that is essentially overlaid on top of our current environment. When you activate a conda environment, you still have access to programs that you installed separately to conda. Therefore, when you run code using conda, your environment is not entirely isolated from your host system. docker, on the other hand, allows you to work with containers. Containers are similar to conda environments in the sense that they contain the software and dependencies that you need to run your code. They are different, though, because containers are entirely isolated from your host system. This can make them more consistent, because users will get the same results regardless of their environment or operating system. 