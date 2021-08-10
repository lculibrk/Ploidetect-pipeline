# Ploidetect-pipeline

  

Data pre-procesing pipeline for https://github.com/lculibrk/Ploidetect

  
## 0. Clone

```
git clone https://github.com/lculibrk/Ploidetect-pipeline
cd ./Ploidetect-pipeline/
```

## 1. Configuration

  

### 1.1: The CONFIG.txt

  

The CONFIG.txt is a `yaml`-formatted file which specifies a number of things required for running the Ploidetect workflow:

  

- The sample names, library names, and paths to the bam files

- The genome version

- Chromosome names

- Sequence type

- Dependency handling

  

This section will break down how to set it up.

  

#### 1.1.1: Specifying your samples

In the CONFIG.txt you will find the following lines:

 
```

bams:
    COLO829-TestA:
        somatic:
            A36971: /projects/POG/POG_data/COLO829-TestA/wgs/biop1_t_A36971/merge/hg19a_bwa-mem-0.7.6a/A36971_2_lanes_dupsFlagged.bam
        normal:
            A36973: /projects/POG/POG_data/COLO829-TestA/wgs/blood1_n_A36973/merge/hg19a_bwa-mem-0.7.6a/A36973_1_lane_dupsFlagged.bam

```

Let's break down this entry further, line-by-line.

  

- `bams:` designates this as the block of samples and bams. **Do not change**

- `COLO829-TestA:` is the **sample name**. For each case (tumor/normal pair), you have some kind of sample name. For example, patient ID could go here.

- `somatic` and `normal` designate the somatic and matched normal libraries. **Do not change**

- `A36971` and `A36973` are the **library names**. Change these to a unique identifier for each sequenced library. They must be unique (don't re-use library names across different samples)

- `/projects/yadda/yadda.bam` is the **path to the bam file** for the somatic and normal libraries. Under `somatic:` you will enter `libraryname: /path/to/the/somatic/bam/file.bam`, and do the same for the normal bam under `normal:`.

  

This sounds a bit annoying, but it's not too bad once you've done it once or twice. If you have a ton of samples we would recommend writing some script to create this automatically and `cat` it into a CONFIG.txt without the `bams` section.

  
  

#### 1.1.2: Specifying the hg version

  

Next there's the genome name section, so change it to match the genome name. So far only hg19 and hg38 are officially supported, but general-genome support is in the works.

  

```

# genome_name should match bams
genome_name:
    "hg19"

```

  

Change the "hg19" in the quotations to whatever you're using.

  

#### 1.1.3: Ploidetect version

This only applies if you're not running the pipeline with `--use-singularity` (to be explained).

  

Specify the version name of Ploidetect to use from the [github repo](https://github.com/lculibrk/Ploidetect) in the `ploidetect_ver` entry. Here we chose v1.0.0 (effectively the paper version), although we've gone up to v1.2.2 so far. If you have a local clone on your filesystem, you can specify that using `ploidetect_local_clone`. If you don't have a clone, leave it blank (i.e. the line should simply read `ploidetect_local_clone: `)

  

```

# ploidetect_ver should be a branch or tag. Overriden by ploidetect_local_clone.
ploidetect_ver: v1.0.0
# Leave ploidetect_local_clone blank or 'None' to download from github
ploidetect_local_clone: /gsc/pipelines/Ploidetect/{ploidetect_ver}

```

  

#### 1.1.4: Conda or singularity?

  

If running the pipeline with `--use-singularity`, set this to 0. Otherwise set it to 1. This option exists because we cannot `conda install` the Ploidetect package, so this flag tells the workflow whether it needs to install the R package itself or not. If we're using singularity, the package is already installed in the container.

  

```
# are we using docker?
install_ploidetect: 0
```

We recommend using singularity to handle the dependencies. This will be explained in section 1.2.
  

#### 1.1.5: hgver-specific files

  

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
        resources/snp_arrays/hg19/SNP_array_positions.txt
    hg38:
        resources/snp_arrays/hg38/SNP_array_positions.txt
```

  

Under `genome`, select your hg-ver and specify the path (relative or absolute) to the genome fasta file. The fasta index (.fai) must also be there. You can create this using `samtools faidx`, for example.

  

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
 - chr2
 - etc

```

  

Once the config has been dealt with, we then handle dependencies.



## 2. Dependencies

### 2.1:  Snakemake


Ploidetect is most easily run using this Snakemake workflow, which take tumor/normal bam files and runs Ploidetect as well as its data preparation steps.

You will need to install `snakemake` to run this workflow. Installation instructions for Snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). You can also run the following command to create a conda environment in the current working directory named `venv` with Snakemake installed

```

make snakemake_env
source venv/bin/activate

```

### 2.2: Package management

Some form of package management is recommended to satisfy the dependencies of the workflow. Currently, Conda and singularity are supported. You only need one to run this workflow. If you followed the above instructions you will already have a conda environment active.

#### 2.2.1 Optional - Install mamba

mamba can be used inplace of conda, and runs much, much faster.

Simply replace 'conda' with 'mamba' in every command after installing. If you installed conda and Snakemake using the Makefile  (section 2.1) Mamba will already be installed and you can skip this step.

```
conda install mamba -c conda-forge
```

#### 2.2.2: Install only Snakemake and let Snakemake handle the dependencies via Conda.
**NOTE: You only need to do ONE OF 2.2.2 or 2.2.3**

If you did not do step 1.1, run the below command:

```
conda install snakemake --yes
```

Then install the dependencies before running the workflow using the following:

```

snakemake --use-conda --conda-create-envs-only -j 1

```
If you installed mamba, you can use the below command instead:
```

snakemake --use-conda --conda-frontend mamba --conda-create-envs-only -j 1

```

Once complete, proceed to section 3

#### 2.2.3: Use singularity

You don't need to do anything beforehand. The container will be automatically pulled for the run. Simply ensure that Singularity is installed on your system.

 
  

## 3. Running

  

We recommend running this workflow with THREADS set to at least 24, to allow parallel processing of the chromosomes as much as possible. We **highly** recommend dry-running Snakemake once after you're done to ensure you have everything set up correctly by appending `-n` to the end of your command. If you have a screen of green/yellow text, you're good. If there's red, then mistakes have been made.

Depending on how you completed section 2.2, you run the workflow using one of three options:

### 3.1. For those who followed 2.2.2 (Conda)

```

snakemake -p --restart-times 2 -j THREADS --use-conda

```

You may encounter an error to the tune of an R lazy-load database being corrupt. `--restart-times 2` should deal with this, but if it doesn't, it should be solvable by just rerunning the command. if that doesn't, open an issue.


### 3.2. For those who followed 2.2.3 (Singularity)

 
```

snakemake -p --restart-times 1 -j THREADS --use-singularity

```

  

If you encounter errors with singularity mode that involve files being not found, first check for typos in the file paths in question. If they're fine, then you probably need to bind the relevant directory. For example if you're running this under a different "root" directory form the data. ie. the data are stored under `/data` and you're running this in `/analysis` or something along those lines. In this situation you need to bind the directories:

  

```

snakemake -p --restart-times 1 -j THREADS --use-singularity --singularity-args "-B /path/to/where/the/bams/are/stored"

```

  

## 4. Optional: cluster submission

  

If you're running on a cluster with `slurm` installed, you can have Snakemake submit the jobs like so:

  
### Conda (2.2.2):
```

mkdir logs_slurm/

snakemake \
 -p \
 --restart-times 2 \
 -j THREADS \
 --use-conda \

```
### Singularity (2.2.3):
```

mkdir logs_slurm/

snakemake \
 -p \
 --restart-times 2 \
 -j THREADS \
 --use-singularity \
 --cluster 'sbatch -t 2-00:00:00 --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{wildcards} -e logs_slurm/{rule}_{wildcards} -p QUEUE_NAME'

```

## 5. Miscellaneous: GSC runner

  

Use GSC bioapps details to automatically find configuration details and run the pipeline for a given ID and biop

  

Eg. run Ploidetect-pipeline for POG965 biop2.

  

```

snakemake -s Snakefile.gsc.smk -pr --config id=POG965 biopsy=biop2

```
