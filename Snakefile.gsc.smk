#! snakemake
"""Run Ploidetect-pipeline with GSC sample setup helpers."""
import os
import sys
from types import SimpleNamespace

sys.path.insert(0, workflow.basedir)
scripts_dir = os.path.join(workflow.basedir, "scripts")
sys.path.append(scripts_dir)
from constants import VERSION as pipeline_ver
from gsc_build_config import CONFIG_BASENAME, get_biopsy_dna_tumour_normal, get_gsc_output_folder
from gsc_build_config import main as build_config

USAGE = """\
Run Ploidetect-pipeline with GSC sample setup helpers.

Use --config options such as:
    id: (required)
    biopsy: eg biop2 - optional - find libraries
    tumour_lib:
    normal_lib:
    gsc_config_filename: (optional) specify output config filename
    project: (optional) specify project

Examples:
    Snakefile.gsc.smk --config id=POG965 biopsy=biop2

    Snakefile.gsc.smk --config id=COLO829-TestA tumour_lib=A36971 normal_lib=A36973 project=POG
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

# Current output_dir value, before loading any defaults.
output_dir = config["output_dir"] if "output_dir" in config.keys() else ""

# Load default config values.
configfile: os.path.join(workflow.basedir, "CONFIG.txt")

output_dir = output_dir if output_dir else get_gsc_output_folder(
    patient_id=config["id"],
    tumour_lib=config["tumour_lib"],
    normal_lib=config["normal_lib"],
    pipeline_ver=pipeline_ver,
    ploidetect_ver=config["ploidetect_ver"],
    project=config["project"] if "project" in config.keys() else None,
)

# Find an output filename for the config we are generating
if "gsc_config_filename" in config.keys():
    output_filename = config["gsc_config_filename"]
else:
    output_filename = os.path.join(output_dir, CONFIG_BASENAME)

# Script building the config if needed
if not os.path.exists(output_filename):
    print(f"Creating {output_filename}")
    args = SimpleNamespace(**config)
    args.pipeline_ver = pipeline_ver
    args.output_file = output_filename
    args.patient_id = args.id
    build_config(args=args)

# Run the standard Snakemake pipeline with the appropriate config
print(f"config: {os.path.abspath(output_filename)}")
config = dict()  # Remove any existing values
configfile: output_filename

# Setting the workdir is less flexible, but should keep logs and paramters organized.
workdir: output_dir

module ploidetect:
    snakefile: "Snakefile"
    config: config

use rule * from ploidetect
