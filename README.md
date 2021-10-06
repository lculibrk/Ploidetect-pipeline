# Ploidetect-pipeline

  

Data pre-procesing pipeline for https://github.com/lculibrk/Ploidetect

  
## 0. Clone

```
git clone https://github.com/lculibrk/Ploidetect-pipeline
cd ./Ploidetect-pipeline/
```

## 1. Configuration

  

### 1.1: The configuration files

  

Under `config/` are three `yaml`-formatted files which specify some default settings, run parameters, and sample information.

  

This section will break down how to set it all up.

  

#### 1.1.1: Specifying your samples

In the `config/samples.yaml` you will find the following lines:

 
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

  

This sounds a bit annoying, but it's not too bad once you've done it once or twice. If you have a ton of samples we would recommend writing some script to create this automatically. In the future we plan to include a script to generate this from a tab-separated file.


#### 1.1.2: Genome version considerations

We're moving on to the `config/parameters.yaml` section

Next we specify the genome name in the config.
  
```
# genome_name should match bams
genome_name: "hg19"
```

Change the `"hg19"` in quotations to whatever you're using.

We also allow you to specify custom cytobands if you are not using hg19 or hg38, or if you want to use a different annotation than those on UCSC. 

```
cyto_path: auto
```

Leave this on "auto" if using hg19 or hg38. If you'd like to use custom cytobands, this must point to a tab-separated file with five columns: chromosome, start, end, band name, and band annotation. Column names are optional. 

#### 1.1.3: Ploidetect version

This only applies if you're not running the pipeline with `--use-singularity` (to be explained).

Specify the version name of Ploidetect to use from the [github repo](https://github.com/lculibrk/Ploidetect) in the `ploidetect_ver` entry. Here we chose v1.3.0 (latest version). If you have a local clone on your filesystem, you can specify that using `ploidetect_local_clone`. If you don't have a clone, leave it blank (i.e. the line should simply read `ploidetect_local_clone: `). You can also specify a commit id, such as the one in the paper. 

  

```

# ploidetect_ver should be a branch or tag.  Overriden by ploidetect_local_clone.
ploidetect_ver: v1.3.0

# Leave ploidetect_local_clone blank or 'None' to download from github
ploidetect_local_clone: 

```

#### 1.1.4: Sequence type
We've recently added support for oxford nanopore data. Specify here whether your data are from short reads (Illumina, BGI) using `short` or from long reads (Oxford nanopore) using `ont`. While we don't foresee issues with using Pacific Biosciences data, it has not been tested. 
```
sequence_type: "short"
```
#### 1.1.5: Supporting custom genomes

Currently, we support hg19 and hg38 without requiring any user-specified files. Genomes hosted on UCSC aside from hg19 and hg38 are compatible, requiring only a few files from the user. If you are using hg19 or hg38, but would like to use your own versions of these files, you may also specify them here. By default the workflow will use packaged files (SNP positions) or UCSC data (fasta and cytoband annotations). 

If you are using a non-hg19/38 genome, you must include these lines in the config:
  
```
array_positions:
    genome_name:
        /path/to/snp/positions.txt
```

This file must be formatted as the files under `resources/snp_arrays/*/SNP_array_positions.txt`, and include the positions of known SNPs in your genome. We used about 200,000 SNP positions to develop and benchmark Ploidetect. Your mileage may vary if you use a different number relative to the length of your genome.

Next, specify a path to your genome fasta file. The directory of the fasta must also contain the fasta index (.fai). You can generate this using `samtools faidx`, for example.

```
genome:
    genome_name:
        /path/to/genome/fasta.fa
```
  
Finally, you must specify the chromosomes to be analyzed. For hg19 and hg38, this defaults to chr1-22 + X. Below is an example of the default hg19 chromosome configuration. The chromosome names must match the .bam that you use. 
  
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
