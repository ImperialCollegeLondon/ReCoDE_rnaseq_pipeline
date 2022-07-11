nextflow pipeline
=================

# What is nextflow?

nextflow is a workflow management system that allows users to orchestrate the running of scripts into a single cohesive pipeline. nextflow pipelines allow reproducibility and can be easily shared and used by others. nextflow is a DSL (domain specific language) with simple methods for combining scripts into fully-fledged pipelines. For detailed information, see [the nextflow site](https://www.nextflow.io/). 

In this document we will discuss features of nextflow that are useful for coordinating our simple RNA-seq pipeline, but these explanations will be brief and will only scratch the surface of what nextflow can do. We will explain the main concepts alongside the building of our simple RNA-seq pipeline. If you want a more in-depth explanation of writing nextflow code, you can read the [nextflow documentation](https://www.nextflow.io/docs/latest/index.html) or look for more detailed guides.

![Example nextflow code](https://www.nextflow.io/img/home-dsl2.png "Example nextflow code, available from https://www.nextflow.io/.")

nextflow can integrate lots of different types of scripts, including those written in bash, Python and R. So, you aren't required to fully learn a new language in order to create nextflow pipelines, you can instead reuse existing scripts and workflows. Our previous pipelines both made use of the same scripts in the `bin/` directory, but they used them in slightly different ways. In this document we will use these same scripts, but we will orchestrate them using nextflow.

When we wrote the simple local pipeline and the parallelised pipeline, we had to modify the scripts in order to make them compatible with the Imperial cluster. If we were to send this pipeline to a researcher at another university, they might have some difficulty adapting it. nextflow is aware of common cluster schedulers, like PBS and slurm, so it can worry about queueing jobs on the cluster for us! We can write code that runs on our local computer, then export it to a computing cluster to do the heavy lifting without changing anything about the pipeline itself. This is possible because we can give nextflow a config file specifying how it can interact with our specific computing cluster. Many institutions have [pre-built config files](https://github.com/nf-core/configs), so you may not even need to write your own!

For our pipeline, we are largely motivated to use nextflow because: i) we don't want to have to worry about setting up our scripts as jobs on the cluster, like we did in `docs/parallelised_pipeline.md`; ii) setting up our pipeline with nextflow will make it easier for others to run it, especially if they use a different computing cluster; iii) nextflow will make sure each processing stage produces the expected outputs, and make it far easier to identify the source of the problem. 

# Environment management

We will also see that nextflow makes it easier to work with environment managers like conda. We have already been working with conda, so hopefully you have picked up the basics. conda lets us create isolated environments that include the packages we need and their dependencies. conda lets us more easily share our code with others, because we can share our code alongside a list of software packages that are required. The easiest way to do this is with an `environment.yml` file, like the one we have created for this exemplar:
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

If we run the code `conda create env -f environment.yml`, like we did at the start of this project, it creates an environment with the name "recode_rnaseq". To create this environment, conda installs all the packages listed as dependencies and makes sure that all of the package versions are compatible with each other. conda knows where these packages are located because we supply the channels `bioconda`, a repository of bioinformatics software packages, and `conda-forge`, another repository for open source packages.

conda has some limitations though. In our case, we want to create a reproducible pipeline that can be run easily by anyone. For instance, we might find that less common packages are not available on conda. Or, we might find that the packages we want are only available on certain operating systems; STAR, for instance, cannot run on Windows. We can avoid some of these limitations with docker, that is also well-integrated with nextflow. 

conda allows us to install a set of packages into an environment that is essentially overlaid on top of our current environment. When you activate a conda environment, you still have access to programs that you installed separately to conda. Therefore, when you run code using conda, your environment is not entirely isolated from your host system. docker, on the other hand, allows you to work with containers. Containers are similar to conda environments in the sense that they contain the software and dependencies that you need to run your code. They are different, though, because containers are entirely isolated from your host system. This can make them more consistent, because users will get the same results regardless of their environment or operating system. 

## Using containers

Previously, we created a single conda environment in which we installed all of the packages we needed to run the pipeline. For the nextflow pipeline, we will instead create an environment for each step of the pipeline. By this, we mean that we will separate each step (QC, trimming, indexing, alignment, counting, multiqc) and use an environment that contains only the package we need. For instance, in the QC step we will create an environment that contains only FastQC and its dependencies. Using different environments for each step is useful as the project gets bigger, so that we don't have to create huge and complex conda environments. 

As we will see when we create the nextflow pipeline, we don't need to set up these environments manually; nextflow will create the environments for us! conda integration is a relatively new feature to nextflow, and it is generally recommended to user docker containers instead. So, we need to make a container for each of the pipeline steps that contains the applications we need. nextflow is helpful in that it will download and activate the containers we specify automatically, but these containers must exist to begin with!

Luckily for us, containers have already been created for many common bioinformatics tools. They are created and hosted by the [BioContainers project](https://biocontainers.pro/). Try visiting their website and finding a container image for the FastQC program. When we write our pipeline we will show you how to specify the container image that contains the relevant packages. First, though, we will discuss how to create your own containers. We have already found a FastQC container image that we could use, but for the sake of learning we will try and create our own FastQC container.

Containers are stored as "images", which are instructions on how to assemble a docker container. We can create our own image that contains FastQC, and use this to build our environment as part of our pipeline. After we have downloaded docker and created a dockerfile (you can find ours in the top-level directory of the exemplar), we can begin specifying our image. To do this, we will be creating an image that is similar to the BioContainers FastQC image. 

The first step in creating the image is to choose a base image. The base image is the docker image you will build your own image on top of. For instance, you may choose to use a minimal linux image to build upon. We will use the base biocontainers image, which uses Ubuntu as its operating system and includes conda within the container. We do this like so:
```
FROM biocontainers/biocontainers:v1.1.0_cv2
```

Below, we use `USER root` to switch to the root user, because we need admin permissions to download and setup FastQC. We then download FastQC and its dependency java, unpack it and allow it to be executed. 

```
USER root

RUN apt-get update && apt-get install -y openjdk-8-jre-headless
RUN mkdir -p /opt/fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip --no-check-certificate -O /opt/fastqc/fastqc_v0.11.9.zip
RUN unzip /opt/fastqc/fastqc_v0.11.9.zip -d /opt/fastqc
RUN rm /opt/fastqc/fastqc_v0.11.9.zip
RUN mv /opt/fastqc/FastQC/* /opt/fastqc
RUN rmdir /opt/fastqc/FastQC && chmod +x /opt/fastqc/fastqc
```

Next, we create a link to the program in `/usr/local/bin/` and add the directory to our path. Now, when we run the `fastqc` command within our container, it will run the program.
```
RUN ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc
ENV PATH /usr/local/bin:$PATH
```

We would prefer that the container be run with a non-root user. This is considered best practice for security: why run FastQC as an admin when we can do it as a standard user? So, we create and switch to a non-root user for when we actually run the container.
```
RUN adduser --system --no-create-home nonroot
USER nonroot
```

To use our docker container, we first have to build it as an image. To do this, we will use the `docker build` command. For this exemplar, I created an account on dockerhub, which allows users to host and share docker images, and stored the image on dockerhub. So, after building the image, I additionally pushed the image to dockerhub, as so:
```
# docker build -t jackgisby/fastqc .
# docker push jackgisby/fastqc
```

This way, you can use my docker image or edit the pipeline to use the image you have created. Usually our next step would be to demonstrate how to access and use the docker image directly, but nextflow will take care of this for us, as we will demonstrate in the following sections. If you want to use docker in your other work, it may be useful to read the [docker documentation](https://docs.docker.com/) or a specific course on docker.

# Running the nextflow pipeline

We can now start building our nextflow pipeline. We will create a module for each of our processing steps in `bin/` and connect them together into a cohesive nextflow pipeline. We will allow users to use either docker or conda to run the pipeline. We will test the pipeline locally before running it on the cluster for the full dataset. Finally, we will look at how we can set up GitHub actions to test the pipeline every time we make a change. 

We will go over the main features of our workflow, but we will not spend much time providing background to nextflow. It may be useful to consult other resources, like the [nextflow documentation](https://www.nextflow.io/docs/latest/index.html), to find out more details or clarify how elements of the pipeline work.

## Building the pipeline processes

The first step in creating our pipeline will be to create a nextflow `process` to run each of the steps in `bin/`. We will create these in the `modules/` directory. For instance, a minimal process for running FastQC on a `.fastq` file might look like this:
```
process FASTQC {
    input:
    path fastq

    output:
    path "*_fastqc.html", emit: fastqc_html
    path "*_fastqc.zip",  emit: fastqc_zip

    script:
    """
    fastqc -o ./ ./$fastq
    """
}
```

The key parts of this process definition include:

 - **input** - The inputs are the files that are going to be used as input to this process. In this case we want to process a `.fastq` file, so we specify that we expect a path as input, and we name this "fastq".

 - **output** - The outputs are the files that we expect our program is going to generate. We don't need to specify all of the program's outputs, just the ones that we want to keep. In this case, we want to capture the files that end in "_fastqc.html" and those that end in "_fastqc.zip". These represent the FastQC report that we will later pass to MultiQC in order to view the QC metrics for our `.fastq` files. We use `emit:` to name the outputs to make accessing them easier later. We use the "*" symbol as a wildcard, meaning this symbol can match anything. This means the 

 - **script** - The script contains one or more commands to be run as part of the process. This is like a bash script, and it can access the variables we defined above in the process. For instance, we use the `fastq` variable to access the location of the `.fastq` file that needs to be processed. This script will run `fastqc` on this `.fastq` file and save the results to the working directory. 

We will discuss later how to pass the inputs and recover the outputs from these processes. You might notice that the actual module, `modules/fastqc.nf` is slightly more complex, demonstrated below:
```
process FASTQC {

    label "long"
    publishDir "$params.outdir/a_fastqc/", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)

    // use a custom container rather than biocontainers
    // container "quay.io/biocontainers/fastqc:0.11.9--0"
    container "jackgisby/fastqc:v1"

    input:
    tuple val(accession), path(fastq)

    output:
    path "*_fastqc.html", emit: fastqc_html
    path "*_fastqc.zip",  emit: fastqc_zip

    script:
    """
    mkdir fastqc_tmp

    $baseDir/bin/fastqc.sh \
        "./" \
        "$fastq" \
        "./fastqc_tmp"
    """
}
```

One major addition is that our process now accesses variables stored within `params`. This variable contains parameters that let us modify how our pipeline runs. We will discuss how we can set and modify these parameters later.

The key additions to the module are as follows:

 - **label** - We can add our processes to different groups using labels. We will make use of labels later when we use them to decide which cluster resources we will allocate to each process.

 - **publishDir** - Even though we have specified expected outputs, nextflow does not necessarily make it easy for us to retrieve our results. nextflow stores all of the input and output files into a directory named `work/`. This is where nextflow does its processing, and it makes it easier for nextflow to isolate each of the processes when they may be running concurrently. If we specify `publishDir`, nextflow will copy the expected output files to a directory of our choosing. We specify the directory as `$params.outdir/a_fastqc/` - i.e. we want the results saved to the `a_fastqc` directory within the folder specified by `$params.outdir`. We set the `mode` variable which controls the method nextflow uses to make the results available.

 - **conda** - If we are using conda to get the FastQC program, then this tells nextflow where it can find FastQC. The variable `params.enable_conda` variable indicates whether we are using conda. If so, we tell nextflow to use a specific version of FastQC available from Bioconda ("bioconda::fastqc=0.11.9"), else we use `null` to indicate nextflow should not download a conda environment.

 - **container** - If we are not using conda, we are going to use containers instead. This command specifies which container nextflow should set up to run the analysis. To do this, we could use the container "quay.io/biocontainers/fastqc:0.11.9--0", which is the URL of the biocontainers FastQC image. But, for this process, we will instead use the docker image we created earlier: "jackgisby/fastqc:v1". This will get version 1 of the FastQC container I published to dockerhub.

 - **input** - We changed the input so that it accepts both the path of the `.fastq` file, but also its sample ID. This might be useful for naming the outputs using the variable `accession`. 

 - **script** - Instead of running FastQC directly, we are now using the script `bin/fastqc.sh` instead. In general I don't recommend doing this, but we will do it this way for consistency with the previous iterations of the pipeline that we developed! Remember that nextflow can also run other types of scripts - in this section we could alternatively run a Python or an R script. Note that we use the `baseDir` variable: this is the base directory of the nextflow pipeline (i.e. the top-level directory of this repository). 

We have created similar processes in the `modules/` directory for each of our pipeline steps. Many of these are very similar, we just run the relevant script in `bin/` while specifying the input and output files that are specific to the program that has been run. 

Next we will move on to specifying the sequence of the analysis, so it might be worthwhile to look at the other processes and try and understand what is going on. You should be able to match up the inputs and outputs to the files we generated for the previous pipelines. 

## Putting the pipeline together

We will now create a script, `workflows/nextflow_pipeline.nf`, to create a cohesive workflow. This will simply join the processes we have created into a single script. With nextflow, we don't explicitly state the order tasks should be run. We simply define how the data flows between each of our processes, and nextflow will decide on the running order and will run steps in parallel where necessary!

The first part of our workflow script simply imports the processes that we generaed previously, for example:
```
include { FASTQC } from "$baseDir/modules/fastqc.nf"
```

Instead of defining a process, we now define a workflow. Workflows are used for orchestrating the flow of data between  multiple processes. A more complex project might even define workflows within other workflows, but our use case is a lot simpler. We define a workflow called `PROCESS_RNASEQ` that will represent our end-to-end pipeline.
```
workflow PROCESS_RNASEQ {
}
```

Within this workflow, we start by getting the input data. nextflow uses variables, called "channels", to describe how data will flow through a workflow. There are two different types of channels. "Value" channels are like regular variables we would define in bash, Python or R. It contains a single value and can be passed to nextflow processes to change how they run. The other type of channel is a "Queue" channel. These channels typically contain multiple inputs that each need to be processed. 

There are multiple ways of defining these channels, but we will use `fromPath`. This allows us to create a channel from a set of input files, like our `.fastq` files. To create our queue, we define a channel called `input_fastq` using `Channel.fromPath` with the `.fastq` files as input. We specify that the channel should look into the `params.fastqc_dir` directory and include any files within it that end with ".fastq.gz" using the "*" wildcard character. This on its own is enough to setup a channel containing all of the `.fastq` files. We do an additional step, using the `.map` command, that adds the name of the sample to the channel. In the process modules, this lets us get both the sample name and path.
```
input_fastq = Channel.fromPath("${params.fastqc_dir}/*.fastq.gz").map {
    tuple( it.name.split('\\.')[0], it )
}
```

In the next step, we use the `FASTQC` as if it was a function by running: `FASTQC(input_fastq)`. nextflow understands that, since there are multiple files in the `input_fastq` channel, it must run the `FASTQC` process for each element of the channel. This is the behaviour we want, because we need to run `fastqc` on each sample individually. 

Note that, just because we call the `FASTQC` process first, doesn't mean this will be the first thing nextflow runs. For instance, the first process launched might be STAR indexing. The first process launched won't be the alignment step, though, because nextflow knows that this step has to be launched after the alignment stage. 

In the pipeline's parameters, which we will discuss shortly, we have included `params.trim` and `params.align` parameters. In the workflow, we use an if statement similar to the following:
```
if (params.align) {
    ALIGN(input)
}
```

Therefore, the users of the pipeline could choose to just do the FastQC stage, or they could skip the trimming stage and do the alignment stage on the raw `.fastq` files. For understanding and testing the pipeline, we will keep these parameters equal to `true`, so that we run the entire pipeline every time.

If `params.trim` is `true`, we run the trimming step on the `input_fastq` channel, just like we did with the `FASTQ` process. When we run a process, like `TRIM`, we can access the files that are emitted by the process using `TRIM.out`. So, for instance, the `TRIM` process emits the channel `trimmed_fastq`, which we can access as follows: `TRIM.out.trimmed_fastq`. Using this, we can control the flow of data between our processes.

Following calling `TRIM`, we define a set of channels, like so:
```
if (params.trim) {

    // trim the raw reads and re-apply QC
    TRIM(input_fastq)

    // define channels
    ch_trimmed_fastqc_html = TRIM.out.trimmed_fastqc_html
    ch_trimmed_fastqc_zip  = TRIM.out.trimmed_fastqc_zip
    ch_trimming_report     = TRIM.out.trimming_report
    ch_fastq_to_align      = TRIM.out.trimmed_fastq

} else {

    // if we don't trim, we must align the raw reads
    ch_fastq_to_align = input_fastq

    // trimming skipped, define channels as empty for input to multiqc
    ch_trimmed_fastqc_html = Channel.empty()
    ch_trimmed_fastqc_zip  = Channel.empty()
    ch_trimming_report     = Channel.empty()
}
```

This just copies the channels defined by `TRIM.out` to a set of new variables. Usually, this would not be necessary. We only do it because of the possibility that `params.trim` is set to `false`. If it was `false`, `TRIM.out` would never be defined. But, we need to send the values of `TRIM.out` to MultiQC later in the pipeline. If we directly sent the channels defined by `TRIM.out` to MultiQC, there would be an error if we did not run the trimming step. To get around this, we redirect the output of `TRIM.out` to a set of new variables, and pass these to MultiQC. Then, if `params.trim` was `false`, we can set these channels to empty so that the channels are defined regardless of whether trimming was performed or not.

Note also that, if trimming was run, we would set `ch_fastq_to_align` to be the trimmed `.fastq` files. These are then used as input to the alignment stage later on. But, if trimming is not run, we instead set this variable to be the raw `.fastq` files defined by `input_fastq`. So, regardless of whether trimming was run, `ch_fastq_to_align` will be a queue channel containing a set of `.fastq` files to be processed!

We next use a similar pattern to run alignment and trimming only if `params.align` is `true`. If it is `false`, we define empty channels to the outputs of these processes that are expected by MultiQC. 

The first step within the if statement is the indexing step:
```
STAR_INDEX(file("$params.genome_fasta"), file("$params.genome_gtf"))
```

To this channel, we use the `file` function to provide two inputs. This function creates a value channel from the path of a single file, which is defined within `params` in this case. The indexing step requires the `.fasta` and `.gtf` files to create the genome index. Alternatively, wee could have assigned each of these value channels to a variable, then passed these variables to `STAR_INDEX`.

We then provide the indexed genome created by this process to `ALIGN`, which creates a `.bam` file that details how the reads were aligned to the genome. The program htseq-count uses the `.bam` alignment file produced by `ALIGN` and the `.gft` annotations file produced by `STAR_INDEX` to create the final `.counts` file that we will use in our downstream analysis.
```
ALIGN(
    ch_fastq_to_align, 
    STAR_INDEX.out.indexed_genome
)

COUNT(
    ALIGN.out.aligned_bam, 
    STAR_INDEX.out.annotation
)
```

The final step is to run MultiQC on the outputs of FastQC, Trim Galore, STAR and htseq-count. At this stage, the nextflow pipeline differs from the previous pipelines. Previously, we ran MultiQC on the FastQC output, then we ran it on the Trim Galore output, then we ran it on the output of STAR and htseq-count. In the nextflow pipeline, we run MultiQC once at the end of the pipeline. This is because users can specify `params.trim = false` and/or `params.align = false`. These outputs are only included in the MultiQC report if they were actually carried out by the pipeline. So, depending on the input parameters, the pipeline can create these same three MultiQC reports. 

As a first step, we recommend trying to run the pipeline on the test dataset on your local machine. This should be easy to accomplish. As long as you have downloaded nextflow and either of conda or docker, you don't even need to have downloaded this repository to run the pipeline! If you run the following code, nextflow will automatically download the repository and run the pipeline for the test dataset using docker!
```
nextflow run ImperialCollegeLondon/ReCoDE_rnaseq_pipeline -profile test,docker
```

You can replace the "docker" profile for the "conda" profile if you prefer conda. Alternatively, if you don't want nextflow to re-download the repository, you can run the following:
```
nextflow main.nf -profile test,docker
```

If you have been using the cluster to run everything so far, you could instead use the above command to run nextflow within a PBS job.

The `main.nf` file located in the top-level directory of the exemplar tells nextflow how to run the pipeline. It simply runs the `PROCESS_RNASEQ` workflow we discussed in this section. Next, we will discuss the config files that are used to configure the pipeline (including defining variables such as `params`) before trying to run the pipeline on the cluster!

## Config files

So far, we have glossed over variables such as `param` that define how our pipeline runs. nextflow pipelines can be configured in a few different ways, including by passing variables via the command line. We will focus on the use of config files to define the behaviour of the pipeline. The main configuration file can be found in the top-level directory by the name of `nextflow.config`. In this file, we define default parameters for all of our `params` and set up the "profiles" like "test" and "docker".

The first section of the file defines all of our parameters. If we were to use the code below, we would define a single parameter, `example_variable`, who has a value of one. 
```
params {
    example_variable = 1
}
```

In our actual config file, we define a few different options that change how our pipeline runs. For instance, it contains the variables `align` and `trim` that determine which sections of the pipeline are run. It also contains variables that let us change the dataset we want to process.

Later in the config, we define `profiles`. The miniature `profiles` section below defines a single profile, "conda", that activates conda in our processes when it is used. So, we would add the command line argument `-profiles conda` to run the pipeline using conda. The file contains another profile for docker, which we used to run the pipeline in the previous section.
```
profiles {
    conda {
        params.enable_conda    = true
    }
}
```

Also in this configuration file is a manifest, which includes metadata relevant to our pipeline. The file also includes the following line:
```
includeConfig 'conf/base.config'
```

This includes all of the configuration options in `conf/base.config` by default. We use the base configuration file, located in `conf`, to define default process settings. We relegate this to another file to demonstrate a basic set of process settings. It is likely that the user will want to overwrite these, because their data may need different resources to the default test dataset we have been using. In this file, we define the default number of cores (`cpus`), amount of memory (`memory`) and running time (`time`) that each process in `modules/` will be allocated by nextflow. This is similar to how we configured our parallelised pipeline on the cluster!

Another important part of the process configuration is the label-based configurations, that start with `withLabel`. This lets us allocate different resources to each of our processes. All of the processes in `modules/` were given a label that corresponds to one of these resource options. To run FastQC, we used the `short` option, trimming uses the `long` option, and there are special configurations for the STAR processes, which require lots of memory. Another benefit of nextflow is that it can automatically re-run our jobs when they fail. This is useful if, for instance, STAR used more memory than expected. If this is the case, nextflow can re-run the job automatically with more memory, preventing a total pipeline failure. We set the default `maxRetries` to 1 and used `errorStrategy` to tell nextflow only to rerun processes for certain types of errors (e.g. reached max runtime, out of memory).

You may have noticed that the label-based configurations use a function called `check_max`. This is a function that makes sure that the resources used do not exceed the maximum values that have been set in `nextflow.config`. This function is defined at the bottom of `nextflow.config`.

In `nextflow.config` we also created a profile called "test". In the previous section, you may have managed to run the test dataset on your local machine using this profile. In the `profiles` section of `nextflow.config`, we load the `conf/test.config` options if the "test" profile is activated. You can see that this configuration overwrites some of the default parameters in order to run the test dataset. For instance, it specifies that we want to use the smaller test dataset. And, it reduces the maximum available resources so that it can run on a standard computer. 

In the next section, we will get the pipeline running on the cluster for the full dataset. To do this, we will use the "imperial_rcs" profile, which configures nextflow to communicate with the PBS scheduler.

## Running the pipeline on the computing cluster

Finally, we are almost read to run the full nextflow pipeline on a proper RNA-seq dataset! In this section, we will assume that you have downloaded the dataset on the cluster, as per the instructions in `docs/parallelised_pipeline.md` or `data/README.md`. Running the pipeline will lead to the generation of results into your home directory, unless you specify otherwise. The intermediate files produced by nextflow (i.e. the `work/` directory) will be stored in the temporary Ephemeral space. We will also assume that you set up the `recode_rnaseq` conda environment for the parallelised pipeline. 

We first need to configure the pipeline on the cluster. The `nextflow.config` file contains a profile called "imperial_rcs" that, if activated, will include the `conf/imperial_rcs.config` configuration. This includes some key configuration options for working with the Imperial cluster. In `params` we set the results folder to a new folder in our home directory. We also set the paths at which nextflow should expect to find the full dataset. The parameters assume you have downloaded this repository to the cluster and stored the data within the `data/` directory. If this is not the case, you will need to modify this configuration or create your own.

We re-define the process options that we originally setup in `base.config`. The aim of this is to target specific resources on the Imperial computing cluster. For instance, at time of writing, the STAR indexing and alignment processes target the "short8" resource. We also additionally specify the `executor`. Setting the executor tells nextflow how it is supposed to run each step of the pipeline. If we set the executor to "pbspro", it will attempt to submit each step as a job to the PBS scheduler. The default executor, local, is designed for running pipelines locally on a machine with multiple cores. We still use it for the "short" and "long" jobs, because these jobs are less computationally intensive and can be run from the original nextflow job. 

We next use an `executor` block to specify the number of jobs that can be sent to the PBS scheduler at a time by nextflow. Currently, you can't run more than 50 jobs at a time on the Imperial cluster. For this pipeline, we set `queueSize` to 7 because there are very few samples we are going to run. 

Finally, we activate singularity. Singularity is similar to docker, in that it can run containers. In fact, singularity can read and run the docker containers we specified earlier. We use singularity on the cluster instead of docker because it is more compatible with the Imperial cluster. For the purposes of running our pipeline on the cluster, we can consider them to be synonymous. We have to specify `autoMounts` and provide the location of our home and temporary spaces on the cluster using the `-B` flag. This makes sure that these folders are accessible within our singularity containers. 

Now that we have our profile set up, we can run nextflow. Our last task is to run nextflow itself as a job. nextflow needs to be run as its own job to coordinate the submission of the pipeline steps to the PBS. For our pipeline, we will also run the processes labelled with "short" or "long" within this coordinator job, purely so that we don't submit too many smaller jobs to the scheduler. These will run within the coordinator job since we have set their executor to "local".

You run the following command to start running the nextflow job on the cluster:
```
qsub nextflow_job.pbs
```

We give the coordinator job 4 cores and 16GB, so that nextflow can be run at the same time as the "local" jobs. The job activates the `recode_rnaseq` conda environment we setup for the parallelised pipeline, as it contains nextflow. We export the following environment variables:
 - **NXF_OPTS** - This ensures Java runs properly on the Imperial cluster.
 - **NXF_TEMP** - This sets nextflow's temporary directory to the cluster's ephemeral space.
 - **NXF_WORK** - This sets the work directory to the cluster's ephemeral space.
 - **NXF_SINGULARITY_CACHEDIR** - This allows singularity to create an image cache in the ephemeral space.

Finally, the pipeline is run to completion. If you have problems running the pipeline, you could try running the test dataset instead to work out what the problem is. Failing this, you could submit the problem as [an issue](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/issues).

## Continous integration

Previously, we tested the pipeline locally before we ran it for the full dataset on the computing cluster. Testing is a great way to ensure that you don't inadvertently break something when you update your code. To take testing another step further, we can set up our tests to run automatically every time we make updates to our code; this practice is known as continous integration. 

One method of setting up continous integration is using GitHub actions. Every time the project code is updated, we update the project on GitHub using a program called git. The main benefit of this process is that git/GitHub track every single change that is made to a project. This allows you to revert things when they go wrong, and you can keep track of all the changes that have been made over time. We can additionally set up a GitHub Actions workflow that will run every time we push changes to GitHub. We have set one such workflow up at `.github/workflows/test_pipeline.yml`.

Including such a file into the `.github/workflows` directory in your repository will be automatically detected by GitHub Actions. It provides instructions on the jobs we want Actions to carry out. For a full explanation of GitHub actions, check out the documentation [provided by GitHub](https://github.com/features/actions). 


In our workflow, we specify `on: [push, pull_request]` to tell Actions that we want the workflow to be run every time we push code to the repository. It will also activate when a pull request is made; pull requests are used to merge new changes into the main project. 

We then specify a single job, run_pipeline, with the name "Run pipeline with test data". Within this job, we designate a strategy in which we create a matrix of parameters describing all the different parameters we want testing. In this case, we have set the job to use ubuntu with the docker profile. We could have included a test on windows by setting `os:` to `[ ubuntu-latest, windows-latest ]`, or, if we wanted to test conda as well, we could have set `profile:` to `[ docker, conda ]`.

This parameter matrix creates a variable called `matrix` that contains the variables of each job. In our case, we only create one job with one set of parameters for our workflow. We indicate the operating system the workflow that will run on using the variable `matrix.os`, and, finally, we include a list of steps that will be run as part of the job. These steps are as follows:

- **Checkout pipeline** - This step uses a pre-made set of instructions available in [another project on GitHub](https://github.com/actions/checkout). It checks out the pipeline and clones it into the current working directory. 

- **Install nextflow** - We download and install nextflow so that we can test the pipeline in a later step. We also set the environment variable `NXF_VER` to use a specific version of the pipeline.

- **Make scripts executable** - By default, the scripts in `bin/` are not executable, meaning they cannot be run within the actions environment. We use `chmod` to ensure the pipeline can run. 

- **Run pipeline** - Finally, we run the pipeline with the test dataset. We set the profile to docker using the `matrix.profile` variable. If nextflow returns an error, for instance if an expected file was not generated, then the workflow will fail. This is useful, because it indicates the code or dependencies have changed in a way that breaks the pipeline.

You can go to [the actions tab on GitHub](https://github.com/ImperialCollegeLondon/ReCoDE_rnaseq_pipeline/actions) to view the runs that have been carried out for each repository push. If you scroll down far enough, you might find runs that failed when the pipeline and tests were incomplete! In this exemplar, we created only a simple workflow for doing a full pipeline test. For more complex pipelines, this may not be sufficient. When pipelines grow, it might be necessary to include separate tests for each segment of the pipeline. This is referred to as unit testing. 

# Limitations

The pipeline is much improved since our first iteration! We can now run the pipeline on the computing cluster, and we can run the processing for each sample in parallel. Furthermore, using nextflow made the pipeline far more robust, because nextflow will handle the submission and tracking of cluster jobs and it will check that all of the expected outputs are created. 

One limitation of using nextflow is that it is slower to setup than simply creating some bash scripts as PBS jobs. But, it probably is a lot faster than manually creating robust pipelines that are both parallelised on the cluster and robust to common errors and failures. 

Another limitation of the pipeline we have created is that it is highly specific to the soybean dataset we are processing. For instance, our pipeline cannot process datasets that use paired end sequencing - this is a very common type of RNA-seq data! Our pipeline also doesn't allow a huge amount of flexibility. We only provide a few basic options in our configuration files. 

We have created an intentionally simple pipeline to demonstrate the concepts of using computing clusters and building nextflow pipelines. If you have some RNA-seq data you want to process with nextflow, we suggest using [the nf-core pipeline](https://github.com/nf-core/rnaseq) which is a highly robust and customisable pipeline. If you wish to gain a better understanding of the pipeline demonstrated in this exemplar, you could try to extend its functionality. You could try to pick a function that the nf-core pipeline has implemented, or choose one of the following:

 - Extend the pipeline to be able to process both single- (like our dataset) and paired-end datasets. Pick a dataset, for instance from the sequence read archive, and try to create an extension of the pipeline to process it. Hint: look into `channel.fromFilePairs`!

 - Add a new tool to the pipeline. For instance, could swap out htseq-count for another tool like featureCounts. 

 - Make the pipeline more customisable. Look up the arguments for one of the programs we have used and try to allow users of the pipeline to change this parameter using a `.config` file.

 - Bonus: The current continous integration tests are very basic. Create a more complex testing procedure that ensures each part of the pipeline is still running correctly. Nextflow's [test](https://github.com/nextflow-io/tests) repository might be helpful for this.

For the creation of this pipeline, we used the program git to keep track of all the changes we made along the way. We then uploaded the project to GitHub so it can be shared with others. If you are interested in learning this too, you could try using git to keep track of the updates you make. The easiest way to do this is to go to the exemplar's GitHub repository and click the "Fork" button. You can now clone your forked repository and begin updating it. Now, when you make changes, you can push them to your forked repository on GitHub.
