# RNA Sequencing and Alignment Pipeline on the WFU DEAC Cluster (DRAFT)

Original workflow: Dr. Danny Arends (https://gist.github.com/DannyArends)

Adapted by: Dr. Sean Anderson (anderss@wfu.edu) and Jessilyn Gao (gaoj20@wfu.edu)


## Introduction

For ease of notation, we will assume that all programs and data reside in a
location called `$RESEARCH`. For this particular example, we can set the value
of that variable like this:

```sh
export RESEARCH=/deac/bio/zhangGrp
```

In other words, when you see `$RESEARCH` or `${RESEARCH}` below, know that it is
exactly equivalent to the path listed there.

To obtain this repo and the scripts clone it into your research path:

```sh
git clone https://github.com/WFU-HPC/zhangGrp-pipeline pipeline
```

and for convenience, you can create a symbolic link to the `scripts` directory

```sh
ln -sf $(realpath ./pipeline/scripts) ${RESEARCH}/scripts
```


## Step 0: Installing the software (ONLY ONCE)

This pipeline uses a variety of open source and free software packages. We can
install these in a centralized location so that all of the scripts have access
to the appropriate packages. **You only need to do the installation once**,
since the software will remain there for subsequent pipeline runs.

The installation script is located at
`${RESEARCH}/scripts/0_installSoftware.sh`, and [follows the original
script](https://gist.github.com/DannyArends/04d87f5590090dfe0dc6b42e5e1bbe15)
closely with some modifications for best practices. You can run this script as
follows:

```sh
bash ${RESEARCH}/scripts/0_installSoftware.sh ${RESEARCH}/software
```

where `bash ${RESEARCH}/scripts/0_installSoftware.sh` is the command to execute
the installation script, and `${RESEARCH}/software` is the location where you
want your software to be installed. The entire installation process should take
around 15 minutes.

In general, you do **not** want to overwrite a previous
installation, so make sure that you are install to a new location, or remove any
previous installation if you want to install to the same location.


## Step 1: Downloading the Genome Files (ONCE PER GENOME)

Now that the software is installed, we are ready to download and build our first
genome. The appropriate script is located at
`${RESEARCH}/scripts/1_buildGenome.sh`, and it can be executed as follows:

```sh
bash ${RESEARCH}/scripts/1_buildGenome.sh \
        ${RESEARCH}/software \
        ftp.ensemblgenomes.org/pub/fungi/release-57/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.chromosome \
        ftp.ensemblgenomes.org/pub/fungi/release-57/gtf/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.57.gtf.gz \
        ftp.ensemblgenomes.org/pub/fungi/release-57/variation/vcf/schizosaccharomyces_pombe/schizosaccharomyces_pombe.vcf.gz \
        ${RESEARCH}/DATA/genome
```

Obviously, this script has a more complicated syntax but we can break it down easily:

```sh
bash ${RESEARCH}/scripts/1_buildGenome.sh \
        "location of installed software" \
        "base URL of chromosome files, excluding the roman numerals and characters after" \
        "complete URL to the .gtf.gz file" \
        "complete URL to the .vcf.gz file" \
        "location of where you want to save the downloaded genome"
```

and note that the `\` symbol merely denotes a line break; i.e. it is all one
single long command broken into several lines for readability.

This command should only take a few seconds, depending on the number of files
that we are downloading. You will get a message when it has completed the
downloads and the files are saved.


### What are we doing?

This script goes to the URL you specify and checks to see which matching files
are available for download; this is why you specify the URL **without** the
roman numerals. This automated approach means that it will work correctly
regardless of the number of genome files that may be available. The GTF and VCF
files are unique, and thus we can just specify the URLs for those files
directly.


## Step 2: Preparing the Genome (ONCE PER GENOME)

Our next step is to prepare the genome from the files that we just downloaded.
The script to do this can be executed as follows:

```sh
bash ${RESEARCH}/scripts/2_prepGenome.sh \
        ${RESEARCH}/software \
        ${RESEARCH}/DATA/genome/Schizosaccharomyces_pombe.ASM294v2.dna.primary_assembly.fa.gz \
        ${RESEARCH}/DATA/genome/Schizosaccharomyces_pombe.ASM294v2.57.gtf.gz \
        ${RESEARCH}/DATA/genome/schizosaccharomyces_pombe.vcf.gz \
        ${RESEARCH}/DATA/genome/STAR
```

which can be broken down in a similar manner,

```sh
bash ${RESEARCH}/scripts/2_prepGenome.sh \
        "location of installed software" \
        "location of primary assembly file" \
        "location of saved .gtf.gz file" \
        "location of saved .vcf.gz file" \
        "location of where you want to save the processed genome"
```

Running this process should take only a few seconds, and will produce some
screen output.


### What are we doing?

This script generates the genome indexes using Samtools, STAR, Tabix, and
Picard. Review the contents of the script to see the individual commands that
are being run.


## Step 3: Leveraging the DEAC Cluster for Batch Processing

The preparation steps detailed above only take a few minutes, but the final
alignment comes at much higher computational expense. Consider that you must
align against each of the samples within your raw data. If it takes two hours to
do one sample, and you have 24 samples -- you are looking at 48 hours of
calculation time if done sequentially!

This is the primary benefit of using the DEAC Cluster -- you can "submit" all of
your calculations simultaneously and finish your entire analysis in a fraction
of the time. Additionally, your calculations will run in *non-interactive* way
that does not require you to be at your computer while they finish. Our goal is
to automate this batch processing of our samples, so that we can submit our jobs
in bulk, and then simply wait for our results to appear within a few hours.

This phase of our work really involves two components:

1. A script that actually runs the alignment pipeline using our datasets, and
2. A script that generates the batch submission commands.

In general, we will not interact directly with 1. and instead use 2. to automate
the process. This can be carried out with the following (very long) command:

```sh
bash ${RESEARCH}/scripts/3a_batchPipeline.sh \
        ${RESEARCH}/software \
        ${RESEARCH}/scripts/3_pipeline.sh \
        ${RESEARCH}/DATA/rna_seq_7.26.2021/raw_data \
        ${RESEARCH}/DATA/genome/STAR \
        ${RESEARCH}/DATA/genome/Schizosaccharomyces_pombe.ASM294v2.dna.primary_assembly.fa.gz \
        ${RESEARCH}/DATA/genome/schizosaccharomyces_pombe.vcf.gz \
        "small" \
        "16" \
        "20G" \
        "00-06:00:00" \
        ${RESEARCH}/DATA/schizosaccharomyces_pombe \
        ./schizosaccharomyces_pombe-BATCH.sh
```

which can be described in detail:

```sh
bash ${RESEARCH}/scripts/3a_batchPipeline.sh \
        "location of installed software" \
        "location of the pipeline script" \
        "location of your raw data (with A1, A2, B1, B2, etc.)" \
        "location of your processed genomes" \
        "location of primary assembly file" \
        "location of saved .vcf.gz file" \
        "Slurm: partition" \
        "Slurm: number of CPU cores" \
        "Slurm: amount of memory with units" \
        "Slurm: requested time in DD-hh:mm:ss format" \
        "desired path to final output directory" \
        "desired path to generated batch submission script"
```

You need to be very careful to provide the correct path to your raw data, which
is the directory that contains the directories like `A1`, `A2`, `B1`, `B2`, etc.

This command does not really *do* anything other than generate a file with the
code to submit all of the jobs. It prints out the final command you should run.
In the example above, the final command could be something like:

```sh
bash schizosaccharomyces_pombe-BATCH.sh
```

**Caution**: running this command will potentially submit dozens of jobs to the
Slurm scheduler. Make sure to review the file first. Please see Appendix 2 below
for a reference of basic Slurm commands for managing your jobs.


### What are we doing?

The pipeline contains all of the commands in the correct order for doing the
alignment. Since there is a specific sequence that needs to be followed, and
each command has a predictable output, we can automate the process from start to
finish. Review the contents of the standalone pipeline script located at
`${RESEARCH}/scripts/3_pipeline.sh` if you wish to understand each part of the
process.


## Step 4: Putting it All Together

A complete example workflow for *Schizosaccharomyces Pombe* can be found at:

```sh
${RESEARCH}/scripts/schizosaccharomyces_pombe-WORKFLOW.sh
```

including software installation and everything.


## Appendix 1: The DEAC Cluster

[The DEAC Cluster](https://hpc.wfu.edu) has almost 100 compute nodes (servers)
with a combined 4000+ CPU cores and 20+ TB of memory. These compute resources
are accessible to all users via the Slurm scheduler which is the "brain" of the
cluster. It allocates resources to user-submitted research jobs and tasks, and
can manage thousands of simultaneous jobs.


## Appendix 2: Basic Slurm Commands

Some basic Slurm commands:

```sh
sinfo -p small      # general state of the nodes (i.e. how busy is the cluster?)

squeue              # view all jobs in the queue
squeue --me         # view all jobs in the queue belonging to you

scancel 4191218     # cancel a specific job with job id 4191218
scancel -u $USER    # cancel ALL jobs belonging to $USER
```


## Disclaimer

This workflow will always have room for improvement, and thus this document may
be constantly revised so check back often. Support for this workflow will be
offered on a "best effort" basis by Dr. Sean Anderson (anderss@wfu.edu).
