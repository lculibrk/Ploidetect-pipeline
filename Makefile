.PHONY: help
.DEFAULT_GOAL := help


VENV_PIP=venv
VENV_CONDA=venv_conda
# Many python versions should work.
PYTHON_VER=3.11


conda_venv:  ## download and install miniconda & mamba. To use 'source venv_conda/bin/activate'
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	bash miniconda.sh -b -u -p $(VENV_CONDA)
	$(VENV_CONDA)/bin/conda install --yes python=$(PYTHON_VER) pip
	$(VENV_CONDA)/bin/pip install -U pip

conda_snakemake:  ## conda install of snakemake into a new '$(VENV_CONDA)'.  After 'source venv/bin/activate'
	# Only snakemake is really needed but other modules are useful for dev.
	$(VENV_CONDA)/bin/conda install --yes -c conda-forge -c bioconda snakemake \
		snakefmt isort black flake8 mypy pylint pydocstyle
	echo "installed conda venv with snakemake in $(VENV_CONDA) activate with:"
	echo "  source $(VENV_CONDA)/bin/activate"

pip_venv:  ## pip install is much faster and smaller, but missing some advanced snakemake features
	# python -m venv $(VENV_PIP)
	# if the local python version is not installed/working use python from conda_venv
	$(VENV_CONDA)/bin/python -m venv $(VENV_PIP)
	$(VENV_PIP)/bin/pip install -U pip

pip_snakemake:  ## python install
	$(VENV_PIP)/bin/pip install -U snakemake snakefmt isort black flake8 mypy pylint pydocstyle

format-code:  ## apply standard formatter like snakefmt and black to scripts.
# python script formatting
	isort --profile black scripts
	black scripts
# Snakefile formatting
	snakefmt Snakefile* scripts

lint: format-code  ## Code formatting and quality checkers
	flake8 --ignore E501 scripts
	mypy --ignore-missing-imports scripts/*.py

lint-snakefiles: format-code  ## snakemake linting suggestions
	snakemake Snakefile.gsc.smk --lint
# must be called individually on each Snakefile, but checks all imported rules
#	snakemake Snakefile --lint

help:  ## show this message and exit
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[32m%-13s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
