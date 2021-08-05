.PHONY: help
.DEFAULT_GOAL := help

help : Makefile
	@sed -n 's/:*##//p' $<

VENV_NAME=venv

snakemake_env: ## conda install of snakemake into a new '$(VENV_NAME)'.  After 'source venv/bin/activate'
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	bash miniconda.sh -b -u -p $(VENV_NAME)
	$(VENV_NAME)/bin/conda install --yes -c conda-forge mamba python=3.8 pip
	# Only snakemake is really needed but other modules are useful for dev.
	# $(VENV_NAME)/bin/mamba install --yes -c conda-forge -c bioconda snakemake
	$(VENV_NAME)/bin/mamba install --yes -c conda-forge -c bioconda snakemake \
		snakefmt jupyter isort black flake8 mypy pylint pydocstyle
	echo "installed conda venv with snakemake in $(VENV_NAME) activate with:"
	echo "  source $(VENV_NAME)/bin/activate"

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
