#! snakemake
"""Run Ploidetect-pipeline with GSC sample setup helpers."""
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, workflow.basedir)
scripts_dir = os.path.join(workflow.basedir, "scripts")
sys.path.append(scripts_dir)
from constants import VERSION as pipeline_ver
from gsc_build_config import (
    CONFIG_BASENAME,
    get_biopsy_dna_tumour_normal,
    get_gsc_output_folder,
)
from gsc_build_config import main as build_config

USAGE = """\
Run Ploidetect-pipeline with GSC sample setup helpers.
Create a config, if it does not yet exist.

Use --config options such as:
    id: (required)
    biopsy: eg 'biop2' (optional) - find libraries tumour_lib & normal_lib by bioapps_api.
    tumour_lib: Required if no biopsy given for bioapps lookup.
    normal_lib:  Required if no biopsy given for bioapps lookup.
    gsc_config_filename: (optional) specify output config filename
    project: (optional) specify project

Examples:
    snakemake -n -s Snakefile.gsc.smk --config id=POG965 biopsy=biop2

    snakemake -n -s Snakefile.gsc.smk --config id=COLO829-TestA tumour_lib=A36971 normal_lib=A36973 project=POG
"""
if "id" in config.keys() and (
    ("biopsy" in config.keys())
    or ("tumour_lib" in config.keys() and "normal_lib" in config.keys())
):
    print("Running GSC build")
else:
    print("No GSC specific config to create.")
    print(USAGE)
    sys.exit()

if "biopsy" in config.keys():
    config["tumour_lib"], config["normal_lib"] = get_biopsy_dna_tumour_normal(
        patient_id=config["id"], biopsy=config["biopsy"]
    )
case = config["id"]
somatic = config["tumour_lib"]
normal = config["normal_lib"]

# Current output_dir value, before loading any defaults.
output_dir = config["output_dir"] if "output_dir" in config.keys() else ""


# Load default values for references / annotations, etc.
configfile: os.path.join(workflow.basedir, "resources/config/default_run_params.yaml")
configfile: os.path.join(workflow.basedir, "resources/config/genome_ref.yaml")


# If no output_dir given, use properties of GSC bioapps_api for standard project output location.
output_dir = (
    output_dir
    if output_dir
    else get_gsc_output_folder(
        patient_id=config["id"],
        tumour_lib=config["tumour_lib"],
        normal_lib=config["normal_lib"],
        pipeline_ver=pipeline_ver,
        ploidetect_ver=config["ploidetect_ver"],
        project=config["project"] if "project" in config.keys() else None,
    )
)

# Find an output filename for the config we are generating
if "gsc_config_filename" in config.keys():
    gsc_config_filename = config["gsc_config_filename"]
else:
    gsc_config_filename = os.path.join(output_dir, CONFIG_BASENAME)

# Script building the config if needed
if not os.path.exists(gsc_config_filename):
    print(f"Creating {gsc_config_filename}")
    args = SimpleNamespace(**config)
    args.pipeline_ver = pipeline_ver
    args.output_file = gsc_config_filename
    args.patient_id = args.id
    build_config(args=args)

# Run the standard Snakemake pipeline with the appropriate config
print(f"config: {os.path.abspath(gsc_config_filename)}")
config = dict()  # Remove any existing values


configfile: gsc_config_filename


# GSC specific - set scratch subfolder as a symlink.
#   A symlink to trans_scratch prevents the project directory from snapshot
#   backups of large temporary files.
scratch = f"{output_dir}/scratch"
if config["temp_dir"] and not os.path.exists(scratch):
    print(f'Creating scratch space:  {config["temp_dir"]}')
    os.makedirs(config["temp_dir"])
if config["temp_dir"] and not os.path.islink(scratch) and not os.path.exists(scratch):
    print(f"Creating symlink: {scratch}")
    os.symlink(config["temp_dir"], scratch)


# Setting the workdir is less flexible, but should keep logs and parameters organized.
workdir: output_dir


container: "docker://lculibrk/ploidetect"


rule complete:
    """Copy number data annotated by genes, has been produced."""
    input:
        f"{output_dir}/cna_genes.bed",


rule annotate_genes:
    """Add Gene info to copynumber change info."""
    input:
        cna=f"{output_dir}/{case}/{somatic}_{normal}/cna_condensed.txt",
        gtf={config["annotation"][config["genome_name"]]},
    output:
        f"{output_dir}/cna_genes.bed",
    log:
        f"{output_dir}/logs/annotate_genes.log",
    conda:
        "conda_configs/r.yaml"
    params:
        scripts_dir=scripts_dir,
    shell:
        "Rscript {params.scripts_dir}/annotate.R -i {input.cna} -a {input.gtf} -o {output} &> {log}"


module ploidetect:
    snakefile:
        "Snakefile"
    config:
        config


use rule * from ploidetect


print(f"Snakefile.gsc.smk target: {rules.complete.input}")
