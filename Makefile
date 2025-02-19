.PHONY: help
.DEFAULT_GOAL := help


VENV_PIP=venv
VENV_CONDA=venv_conda
# Many python versions should work.
PYTHON_VER=3.11
# Some kind of Pulp error with snakemake 8.27
SNAKEMAKE_VER='snakemake<8.2'

conda_venv:  ## download and install miniconda & mamba. To use 'source venv_conda/bin/activate'
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	bash miniconda.sh -b -u -p $(VENV_CONDA)
	$(VENV_CONDA)/bin/conda install --yes python=$(PYTHON_VER) pip
	$(VENV_CONDA)/bin/pip install -U pip

conda_snakemake:  ## conda install of snakemake into a new '$(VENV_CONDA)'.  After 'source venv/bin/activate'
	# Only snakemake is really needed but other modules are useful for dev.
	$(VENV_CONDA)/bin/conda install --yes -c conda-forge -c bioconda $(SNAKEMAKE_VER) \
		snakefmt isort black flake8 mypy pylint pydocstyle samtools
	echo "installed conda venv with snakemake in $(VENV_CONDA) activate with:"
	echo "  source $(VENV_CONDA)/bin/activate"

pip_venv:  ## pip install is much faster and smaller, but missing some advanced snakemake features
	python -m venv $(VENV_PIP)
	# if the local python version is not installed/working use python from conda_venv
	# $(VENV_CONDA)/bin/python -m venv $(VENV_PIP)
	$(VENV_PIP)/bin/pip install -U pip

pip_snakemake:  ## python install
	$(VENV_PIP)/bin/pip install -U pip
	$(VENV_PIP)/bin/pip install -U $(SNAKEMAKE_VER) snakefmt

format-code:  ## apply standard formatter like snakefmt and black to scripts.
	# python script formatting
	$(VENV_PIP)/bin/pip install -U snakefmt isort black flake8 mypy pylint pydocstyle
	$(VENV_PIP)/bin/isort --profile black scripts
	$(VENV_PIP)/bin/black scripts *.py
	# Snakefile formatting
	$(VENV_PIP)/bin/snakefmt Snakefile* *.smk

lint: format-code  ## Code formatting and quality checkers
	$(VENV_PIP)/bin/flake8 --ignore E501 scripts
	$(VENV_PIP)/bin/mypy --ignore-missing-imports scripts/*.py

lint-snakefiles: format-code  ## snakemake linting suggestions
	snakemake Snakefile.gsc.smk --lint
# must be called individually on each Snakefile, but checks all imported rules
#	snakemake Snakefile --lint

help:  ## show this message and exit
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[32m%-13s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
