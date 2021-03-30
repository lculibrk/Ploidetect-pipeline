# Ploidetect-pipeline

Data pre-procesing pipeline for https://github.com/lculibrk/Ploidetect

## Installation

### Install snakemake

Check for a snakemake installation: `which snakemake`

If not found, suggest installing a conda enviroment for snakemake.

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
This option may take longer and be less stable

#### Option2: Install only snakemake

```
mamba install snakemake --yes
```
or
```
mamba install snakemake --yes
```

Now either install for the `--use-conda` or `--use-singularity`/docker options

##### .snakemake conda

```
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only -j 1
```

##### singularity/docker

```
snakemake \
 --use-singularity \
 --singularity-args "-B /path/to/bams,data/,output/" \
 --cluster 'sbatch -t 2-00:00:00 --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{wildcards} -e logs_slurm/{rule}_{wildcards} -p QUEUE'
```
