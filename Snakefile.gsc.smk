import os
import sys

include: os.path.join(workflow.basedir, "Snakefile")

pipeline_ver = __version__  # imported from Snakefile
ploidetect_ver = config["ploidetect_github_version"]

# Somehow determine we want to try autogenerating a config
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
        create_config: (optional) specify output config filename
        project: (optional) specify project
    """)
    sys.exit()

scripts_dir = os.path.join(workflow.basedir, "scripts")
sys.path.append(scripts_dir)
from gsc_build_config import CONFIG_BASENAME, get_biopsy_dna_tumour_normal, get_gsc_output_folder

if "biopsy" in config.keys():
    tumour_lib, normal_lib = get_biopsy_dna_tumour_normal(patient_id=config["id"], biopsy=config["biopsy"])
if "tumour_lib" in config.keys():
    tumour_lib = config["tumour_lib"]
if "normal_lib" in config.keys():
    normal_lib = config["normal_lib"]

# Find an output filename for the config we are generating
if "create_config" in config.keys():
    output_filename = config["create_config"]
else:
    output_filename = get_gsc_output_folder(
        patient_id=config["id"],
        tumour_lib=tumour_lib,
        normal_lib=normal_lib,
        pipeline_ver=pipeline_ver,
        ploidetect_ver=ploidetect_ver,
        project=config["project"] if "project" in config.keys() else None,
    )
    output_filename = os.path.join(output_filename, CONFIG_BASENAME)

print(f"config: {output_filename}")

rule gsc_config:
    params:
        ver=f'--pipeline-ver {pipeline_ver} --ploidetect-ver {ploidetect_ver}',
        id=config["id"],
        libs=f'-b {config["biopsy"]}' if "biopsy" in config.keys() else f'-t {config["tumour_lib"]} -n {config["normal_lib"]}',
        project=f' --project {config["project"]}' if "project" in config.keys() else ""
    output:
        output_filename
    shell:
        "{scripts_dir}/gsc_build_config.py {params.ver} -i {params.id} {params.libs}{params.project} -o {output_filename}"
