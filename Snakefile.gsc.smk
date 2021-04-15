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
from gsc_build_config import get_biopsy_dna_tumour_normal, get_gsc_output_folder, get_ploidetect_temp_folder

if "biopsy" in config.keys():
    tumour_lib, normal_lib = get_biopsy_dna_tumour_normal(patient_id=config["id"], biopsy=config["biopsy"])
if "tumour_lib" in config.keys():
    tumour_lib = config["tumour_lib"]
if "normal_lib" in config.keys():
    normal_lib = config["normal_lib"]



output_dir = get_gsc_output_folder(config["id"], tumour_lib, normal_lib, pipeline_ver, ploidetect_ver, project=config["project"] if "project" in config.keys() else None)
temp_dir = get_ploidetect_temp_folder(config["id"], tumour_lib, normal_lib, pipeline_ver, ploidetect_ver, project=config["project"] if "project" in config.keys() else None)


ruleorder: preseg > gsc_preseg


rule all_gsc:
    input:
        expand("{output_dir}/cna.txt", output_dir=output_dir)

rule gsc_preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        expand("{temp_dir}/merged.bed", temp_dir=temp_dir)
    output:
        "{output_dir}/segmented.RDS"
    resources: cpus=24, mem_mb=189600
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "Rscript {scripts_dir}/prep_ploidetect2.R -i {input} -o {output}"

