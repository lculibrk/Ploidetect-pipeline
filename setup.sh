wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3/
$HOME/miniconda3/bin/conda init
source $HOME/.bashrc
conda install -c conda-forge mamba
mamba install -c conda-forge -c bioconda snakemake
snakemake --use-conda --conda-create-envs-only -j 1

