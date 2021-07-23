# Ploidetect-pipeline

Data pre-procesing pipeline for https://github.com/lculibrk/Ploidetect

## Installation

### First step

`git clone https://github.com/lculibrk/Ploidetect-pipeline`

`cd Ploidetect-pipeline/`

### Dependencies

Ploidetect is most easily run using this Snakemake workflow, which takes care of processing the data from bam to input, and running the R package itself.

You will need to install `snakemake` to run the Ploidetect workflow. From here, some form of package management is recommended to take care of the dependencies.

For dependencies, you have two options: `conda` and `singularity` (containerized). One of the two must be installed.

### Install snakemake

Check for a snakemake installation: `which snakemake`

If not found, you can use the Makefile to create a conda environment with Snakemake installed:
```
make snakemake_env
```

You can also use one of the below options. If you've installed Snakemake, you can skip to the **Configuration** step.

If conda is not installed on the current system, install miniconda:

https://docs.conda.io/en/latest/miniconda.html

#### Optional - Install mamba

mamba can be used inplace of conda, and runs much, much faster.
Simply replace 'conda' with 'mamba' in every command.

```
conda install mamba -c conda-forge
```

#### Option1: Create a single conda environment with everything

```
mamba env create -p env -f gsc_pipeline/conda_configs/ploidetect_all.yaml
```
or conda equivalent command.  -p option can be any path.
```
conda env create -p <path/to/install> -f conda_configs/ploidetect_all.yaml
```
This option may take longer and be less stable. 

Simply activate the environment:

```
conda activate /path/to/created/env
```

Once complete, proceed to the configuration or (if complete), the running section

#### Option2: Install only Snakemake and let Snakemake handle the dependencies.

```
conda install snakemake --yes
```
or
```
mamba install snakemake --yes
```

##### Install conda envs

```
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only -j 1
```
Once complete, proceed to the configuration or (if complete), the running section

## Configuration

### The CONFIG.txt

The CONFIG.txt is a `yaml`-formatted file which specifies a number of things required for running the Ploidetect workflow:

 - The sample names, library names, and paths to the bam files
 - The genome version
 - Chromosome names
 - Sequence type
 - Dependency handling

This section will break down how to set it up.

#### Specifying your samples

In the CONFIG.txt you will find the following lines:

```
bams:
    COLO829-TestA:
        somatic:
            A36971: /projects/POG/POG_data/COLO829-TestA/wgs/biop1_t_A36971/merge/hg19a_bwa-mem-0.7.6a/A36971_2_lanes_dupsFlagged.bam
        normal:
            A36973: /projects/POG/POG_data/COLO829-TestA/wgs/blood1_n_A36973/merge/hg19a_bwa-mem-0.7.6a/A36973_1_lane_dupsFlagged.bam
```
Let's break down this entry further.

 - `bams:` designates this as the block of samples and bams. **Do not change**
 - `COLO829-TestA:` is the **sample name**. For each case (tumor/normal pair), you have some kind of sample name. For example, patient ID could go here.
 - `somatic` and `normal` designate the somatic and matched normal libraries. **Do not change**
 - `A36971` and `A36973` are the **library names**. Change these to a unique identifier for each sequenced library. They must be unique (don't re-use library names across different samples)
 - `/projects/yadda/yadda.bam` is the **path to the bam file** for the somatic and normal libraries. Under `somatic:` you will enter `libraryname: /path/to/the/somatic/bam/file.bam`, and do the same for the normal bam under `normal:`.

This sounds a bit annoying, but it's not too bad once you've done it once or twice. If you have a ton of samples we would recommend writing some python script to create this automatically and `cat` it into a CONFIG.txt without the `bams` section. 


#### Specifying the hg version

Next there's the genome name section, so change it to match the genome name. So far only hg19 and hg38 are officially supported, but general-genome support is in the works.

```
genome_name:
    "hg19"
```

Change the "hg19" in the quotations to whatever you're using. 

#### Versions

This only applies if you're not running the pipeline with `--use-singularity` (to be explained). 

Specify the version name of Ploidetect to use from the [github repo](https://github.com/lculibrk/Ploidetect) in the `ploidetect_ver` entry. Here we chose v1.0.0, although we've gone up to v1.2.2 so far. If you have a local clone on your filesystem, you can specify that using `ploidetect_local_clone`. If you don't have a clone, leave it blank (i.e. the line should simply read `ploidetect_local_clone: `)

```
# ploidetect_ver should be a branch or tag.  Overriden by ploidetect_local_clone.
ploidetect_ver: v1.0.0
# Leave ploidetect_local_clone blank or 'None' to download from github
ploidetect_local_clone: /gsc/pipelines/Ploidetect/{ploidetect_ver}
```

#### Conda or docker?

If running the pipeline with `--use-singularity`, set this to 0. Otherwise set it to 1. This option exists because we cannot `conda install` the Ploidetect package. If we're using containers, the package is already installed in the container.

```
# are we using docker?
install_ploidetect: 0
```

#### hgver-specific files

The next lines are file dependencies that are needed for running the workflow. Currently the `genome` section is probably the only one you have to touch. `annotation` is used for the BCGSC (our institution)-specific module in `Snakefile.gsc.smk`, and the `array_positions` files are included with the repository.

```
genome:
    hg19:
        /gsc/resources/Homo_sapiens_genomes/hg19a/genome/fasta/hg19a.fa
    hg38:
        /gsc/resources/Homo_sapiens_genomes/hg38_no_alt/genome/fasta/hg38_no_alt.fa
annotation:
    hg19:
        /gsc/resources/annotation/Homo_sapiens.GRCh37.87.gtf
    hg38:
        /gsc/resources/annotation/Homo_sapiens.GRCh38.100.gtf
array_positions:
    hg19:
        "resources/snp_arrays/hg19/SNP_array_positions.txt"
    hg38:
        "resources/snp_arrays/hg38/SNP_array_positions.txt"
```

Under `genome`, select your hg-ver and specify the path to the genome fasta file. The fasta index (.fai) must also be there. You can create this using `samtools faidx`, for example. 

The final section deals with chromosome names:

```
chromosomes:
 - 1
 - 2
 - 3
 - 4
 - 5
 - 6
 - 7
 - 8
 - 9
 - 10
 - 11
 - 12
 - 13
 - 14
 - 15
 - 16
 - 17
 - 18
 - 19
 - 20
 - 21
 - 22
 - X
```

This is fairly self-explanatory. If your chosen bams use "chr" prefixes, add them to these lines. e.g.

```
chromosomes:
 - chr1
 - etc
```

Once the config has been dealt with, you can run the workflow!

Now either install for the `--use-conda` or `--use-singularity`/docker options

## Running

We recommend running this workflow with THREADS set to at least 24, to allow parallel processing of the chromosomes as much as possible. We **highly** recommend dry-running Snakemake once to ensure you have everything set up correctly by appending `-n` to the end of your command. If you have a screen of green text, you're good. If there's red, then mistakes have been made.

### Running with conda (option1, with everything in one repo)
```
snakemake -p --restart-times 1 -j THREADS
```

### Running with conda (option2, letting Snakemake handle the dependencies)

```
snakemake -p --restart-times 1 -j THREADS --use-conda
```

### Running with singularity

```
snakemake -p --restart-times 1 -j THREADS --use-singularity
```

If you encounter errors with singularity mode that involve files being not found, first check for typos in the file paths in question. If they're fine, then you probably need to bind the directory. For example if you're running this under a different "root" directory form the data. ie. the data are stored under `/data` and you're running this in `/analysis` or something along those lines. In this situation you need to bind the directories:

```
snakemake -p --restart-times 1 -j THREADS --use-singularity --singularity-args "-B /path/to/where/the/bams/are/stored"
```

### Cluster submission

If you're running on a cluster with `slurm` installed, you can have Snakemake submit the jobs like so:

```
mkdir logs_slurm/
snakemake -p --restart-times 1 -j THREADS --use-singularity \
 --cluster 'sbatch -t 2-00:00:00 --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{wildcards} -e logs_slurm/{rule}_{wildcards} -p QUEUE_NAME'
```

##### singularity/docker with slurm cluster submission

```
mkdir  logs_slurm data output
snakemake \
 --use-singularity \
 --singularity-args "-B /path/to/bams,data/,output/" \
 --cluster 'sbatch -t 2-00:00:00 --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{wildcards} -e logs_slurm/{rule}_{wildcards} -p QUEUE'
```

## GSC runner

Use GSC bioapps details to automatically find configuration details and run the pipeline for a given ID and biop

Eg. run Ploidetect-pipeline for POG965 biop2.

```
snakemake -s Snakefile.gsc.smk -pr --config id=POG965 biopsy=biop2
```
