import os
import sys
from types import SimpleNamespace

sys.path.insert(0, workflow.basedir)
scripts_dir = os.path.join(workflow.basedir, "scripts")
sys.path.append(scripts_dir)
from constants import VERSION as pipeline_ver
from gsc_build_config import CONFIG_BASENAME, get_biopsy_dna_tumour_normal, get_gsc_output_folder
from gsc_build_config import main as build_config

# start with defaults
configfile: os.path.join(workflow.basedir, "CONFIG.txt")
module ploid:
    snakefile: "Snakefile"
    config: config

ploidetect_ver = config["ploidetect_github_version"]

# Check if a GSC specific config exists or must be built
if any((gsc_key in config.keys() for gsc_key in ["biopsy", "tumour_lib", "normal_lib"])):
    print("Running GSC build")
else:
    print("No GSC specific config to create.")
    print("""
    Use --config options such as:
        id: (required)
        biopsy: eg biop2 - optional - find libraries
        tumour_lib:
        normal_lib:
        gsc_config_filename: (optional) specify output config filename
        project: (optional) specify project
    """)
    sys.exit()

if "biopsy" in config.keys():
    config["tumour_lib"], config["normal_lib"] = get_biopsy_dna_tumour_normal(patient_id=config["id"], biopsy=config["biopsy"])

# Find an output filename for the config we are generating
if "gsc_config_filename" in config.keys():
    output_filename = config["gsc_config_filename"]
else:
    output_dir = get_gsc_output_folder(
        patient_id=config["id"],
        tumour_lib=config["tumour_lib"],
        normal_lib=config["normal_lib"],
        pipeline_ver=pipeline_ver,
        ploidetect_ver=ploidetect_ver,
        project=config["project"] if "project" in config.keys() else None,
    )
    output_filename = os.path.join(output_dir, CONFIG_BASENAME)

# Script building the config if needed
if not os.path.exists(output_filename):
    print(f"Creating {output_filename}")
    args = SimpleNamespace(**config)
    args.pipeline_ver = pipeline_ver
    args.ploidetect_ver = ploidetect_ver
    args.output_file = output_filename
    args.patient_id = args.id
    build_config(args=args)

print(f"config: {output_filename}")

# Run the standard Snakemake pipeline with the appropriate config
configfile: output_filename
module ploidetect:
    snakefile: "Snakefile"
    config: config

use rule * from ploidetect
