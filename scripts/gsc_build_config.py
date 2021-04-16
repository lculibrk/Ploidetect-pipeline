#! /usr/bin/env python
"""Build GSC default values and config file."""
import argparse
import glob
import logging
import os
import sys
from datetime import datetime
from os.path import abspath, dirname, exists, join, realpath

from ProjectInfo import BioappsApi
from ruamel_yaml import YAML

__version__ = "0.0.1"
API = BioappsApi()
CONFIG_BASENAME = "Ploidetect-pipeline.yaml"
GENOME_DATA = """\
# Reference data.  Selected by 'genome_name' value.
genome:
    hg19:
        /gsc/resources/Homo_sapiens_genomes/hg19a/genome/fasta/hg19a.fa
    hg38:
        /gsc/resources/Homo_sapiens_genomes/hg38_no_alt/genome/fasta/hg38_no_alt.fa
array_positions:
    hg19:
        "resources/snp_arrays/hg19/SNP_array_positions.txt"
    hg38:
        "resources/snp_arrays/hg38/SNP_array_positions.txt"
chromosomes:
 - 1
 - 2
 - 3
 - 4
 - 5
 - 6
 - 7
 - 8
 - 9
 - 10
 - 11
 - 12
 - 13
 - 14
 - 15
 - 16
 - 17
 - 18
 - 19
 - 20
 - 21
 - 22
 - X
"""


def get_biopsy_dna_tumour_normal(patient_id, biopsy="biop1"):
    """Return DNA (tumour, normal) for the patient biopsy.

    Example:
        >>> get_biopsy_dna_tumour_normal("POG965", biopsy="biop2")
        ('P02866', 'P02590')
    """
    libs = API.get_lib_names(patient_id)
    normals = [
        lib for lib in libs if API.is_normal_library(lib) and API.is_dna_library(lib)
    ]
    tumour_libs = [
        lib for lib in libs if API.get_biopsy(lib) == biopsy and API.is_dna_library(lib)
    ]
    assert len(normals) == 1
    assert len(tumour_libs) == 1
    return (tumour_libs[0], normals[0])


def get_project_info(patient_id):
    """Return a single project_info dict for the patient."""
    proj_info = API.get_biofx_project_info(patient_id, tumour_char_only=True)
    if len(proj_info) > 1:
        logging.error(
            f"Assuming project: {proj_info[0]['name']}"
            f" - ignoring {[pi['name'] for pi in proj_info[1:]]}"
        )
    return proj_info[0]


def get_somatic_folder(patient_id, tumour_lib, normal_lib, project=None):
    """Return the 'somatic' paired analysis folder.

    Example:
        >>> tumour,normal = get_biopsy_dna_tumour_normal("POG965", biopsy="biop2")
        >>> get_somatic_folder("POG965", tumour, normal)
        '/projects/POG/POG_data/POG965/wgs/biop2_t_P02866_blood1_n_P02590'
    """
    project_path = (
        API.get_biofx_project_path(project)
        if project
        else get_project_info(patient_id)["path"]
    )
    glob_str = f"{project_path}/{patient_id}/*/*_t_{tumour_lib}*_n_{normal_lib}"
    somatic_folders = glob.glob(glob_str)
    if not somatic_folders:
        raise ValueError(f"Failed to find somatic folder.  Checked '{glob_str}'")
    elif len(somatic_folders) > 1:
        raise ValueError(
            f"Failed to find a unique somatic folder.  Found {somatic_folders}"
        )
    return somatic_folders[0]


def get_ploidetect_temp_folder(
    patient_id, tumour_lib, normal_lib, pipeline_ver, ploidetect_ver, project=None
):
    """GSC trans_scratch temp folder.

    Example:
        >>> get_ploidetect_temp_folder("COLO829-TestA", "A36971", "A36973", "0.0.1", "v1.0.0", "POG")
        '/projects/trans_scratch/validations/Ploidetect/POG/COLO829-TestA/Ploidetect-pipeline-0.0.1/Ploidetect-v1.0.0/A36971_A36973'
    """
    TEMP_ROOT = "/projects/trans_scratch/validations/Ploidetect"
    if not project:
        project = get_project_info(patient_id)["name"]
    return join(
        TEMP_ROOT,
        project,
        patient_id,
        "Ploidetect-pipeline-" + pipeline_ver,
        "Ploidetect-" + ploidetect_ver,
        f"{tumour_lib}_{normal_lib}",
    )


def get_gsc_output_folder(
    patient_id,
    tumour_lib,
    normal_lib,
    pipeline_ver,
    ploidetect_ver,
    project=None,
):
    """GSC final project output folder.

    Example:
        >>> get_gsc_output_folder("POG965", "P02866", "P02590", "0.0.1", "v1.0.0")
        '/projects/POG/POG_data/POG965/wgs/biop2_t_P02866_blood1_n_P02590/Ploidetect/Ploidetect-pipeline-0.0.1/Ploidetect-v1.0.0'
    """
    somatic_pair_folder = get_somatic_folder(
        patient_id, tumour_lib, normal_lib, project=project
    )
    return join(
        somatic_pair_folder,
        "Ploidetect",
        "Ploidetect-pipeline-" + pipeline_ver,
        "Ploidetect-" + ploidetect_ver,
    )


def get_bam(library, genome_reference=None):
    """Find a bam for library.  Restricted genome_refernce if given.

    Returns (library, genome_reference)
    """
    bams_all = API.get_bam_info_from_library(library)
    if genome_reference:
        bams_all = [bd for bd in bams_all if bd["genome_reference"] == genome_reference]
    bams = []
    for bam in bams_all:
        bam_fns = glob.glob(join(bam["data_path"], "*.bam"))
        if not bam_fns:  # GERO-114 - COLO829 - no bam in path
            logging.error(f"No bam found in path: {bam['data_path']}")
        else:
            bams.append(bam)
    for bam in bams[1:]:
        warn = f"Ignoring tumour bam {bam['data_path']}"
        logging.warning(warn)

    bam_fns = glob.glob(join(bams[0]["data_path"], "*.bam"))
    assert len(bam_fns) == 1, f"Bam finding error for: '{library}'"
    return (bam_fns[0], bams[0]["genome_reference"])


def genome_reference2genome_name(gsc_genome_reference):
    return "hg19" if gsc_genome_reference == "hg19a" else gsc_genome_reference


def build_config(
    patient_id,
    biopsy=None,
    tumour_lib=None,
    normal_lib=None,
    pipeline_ver="undefined",
    ploidetect_ver="undefined",
    use_docker=False,
    project=None,
    **kwargs,
):
    """Build a GSC config for running Ploidetect.

    Checks GSC bioapps api for relevant details of libraries, bams and output paths.

    Example:
        >>> cfg = build_config("POG965", biopsy="biop2")
        >>> YAML().dump(cfg, sys.stdout)  # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
        # Created by: ...
        id: POG965
        libs:
          tumour: P02866
          normal: P02590
        bams:
          somatic: /projects/analysis/analysis30/P02866/merge32550_bwa-mem-0.7.6a-sb/150bp/hg19a/P02866_3_lanes_dupsFlagged.bam
          normal: /projects/analysis/analysis30/P02590/HCW32CCXY_8/P02590/150nt/hg19a/bwa-mem-0.7.6a-sb/P02590_1_lane_dupsFlagged.bam
        genome_name: hg19
        output_dir: /projects/POG/POG_data/POG965/wgs/biop2_t_P02866_blood1_n_P02590/Ploidetect/Ploidetect-pipeline-undefined/Ploidetect-undefined
        temp_dir: /projects/trans_scratch/validations/Ploidetect/POG/POG965/Ploidetect-pipeline-undefined/Ploidetect-undefined/P02866_P02590
        # ploidetect_github_version should be a branch or tag.  Overriden by ploidetect_local_clone.
        ploidetect_github_version: undefined
        # Leave ploidetect_local_clone blank or 'None' to download from github
        ploidetect_local_clone: /gsc/pipelines/Ploidetect/undefined
        use-docker: 0
        # Reference data.  Selected by 'genome_name' value.
        genome:
          hg19: /gsc/resources/Homo_sapiens_genomes/hg19a/genome/fasta/hg19a.fa
          hg38: /gsc/resources/Homo_sapiens_genomes/hg38_no_alt/genome/fasta/hg38_no_alt.fa
        array_positions:
          hg19: resources/snp_arrays/hg19/SNP_array_positions.txt
          hg38: resources/snp_arrays/hg38/SNP_array_positions.txt
        chromosomes:
        - 1
        - 2
        - 3
        - 4
        - 5
        - 6
        - 7
        - 8
        - 9
        - 10
        - 11
        - 12
        - 13
        - 14
        - 15
        - 16
        - 17
        - 18
        - 19
        - 20
        - 21
        - 22
        - X
    """
    TAB = "  "
    yaml_lines = []

    if not biopsy and not (tumour_lib and normal_lib):
        raise ValueError(
            "Either a biopsy or tumour_lib and normal_lib must be supplied."
        )
    elif biopsy:
        tumour_lib, normal_lib = get_biopsy_dna_tumour_normal(patient_id, biopsy)

    prog_str = f'Created by: {realpath(abspath(__file__))} at {datetime.now().strftime("%Y%m%d %H:%M:%S")}'
    yaml_lines.append(f"# {prog_str}")

    # Find bams
    yaml_lines.append("bams:")
    yaml_lines.append(f"{TAB}{patient_id}:")
    tumour_bam_fn, genome_name = get_bam(tumour_lib)
    normal_bam_fn, normal_genome_name = get_bam(normal_lib)

    yaml_lines.append(f"{TAB}{TAB}somatic:")
    yaml_lines.append(f"{TAB}{TAB}{TAB}{tumour_lib}: {tumour_bam_fn}")
    yaml_lines.append(f"{TAB}{TAB}normal:")
    yaml_lines.append(f"{TAB}{TAB}{TAB}{normal_lib}: {normal_bam_fn}")

    assert genome_name == normal_genome_name
    yaml_lines.append(f"genome_name: {genome_reference2genome_name(genome_name)}")

    # output_paths
    yaml_lines.append(
        f"output_dir: {get_gsc_output_folder(patient_id, tumour_lib, normal_lib, pipeline_ver, ploidetect_ver, project=project)}"
    )
    yaml_lines.append(
        f"temp_dir: {get_ploidetect_temp_folder(patient_id, tumour_lib, normal_lib, pipeline_ver, ploidetect_ver, project=project)}"
    )

    # Ploidetect installation and versions.
    yaml_lines.append(
        "# ploidetect_github_version should be a branch or tag.  Overriden by ploidetect_local_clone."
    )
    yaml_lines.append(f"ploidetect_github_version: {ploidetect_ver}")
    yaml_lines.append(
        "# Leave ploidetect_local_clone blank or 'None' to download from github"
    )
    yaml_lines.append(
        f"ploidetect_local_clone: /gsc/pipelines/Ploidetect/{ploidetect_ver}"
    )
    yaml_lines.append(f"use-docker: {1 if bool(use_docker) else 0}")

    # Genomic Reference data
    yaml_lines.append(GENOME_DATA)

    yaml = YAML()
    config = yaml.load("\n".join(yaml_lines))

    return config


def parse_args():
    """Parse command line arguments with argparse."""
    parser = argparse.ArgumentParser(
        description="Build standard config for GSC patient.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--pipeline-ver", required=True, help="Ploidetect-pipeline Version"
    )
    parser.add_argument("--ploidetect-ver", required=True, help="Ploidetect Version")
    parser.add_argument(
        "-i",
        "--patient-id",
        required=True,
        help="GSC bioapps patient_id.  eg POG965",
    )
    parser.add_argument("-b", "--biopsy", help="Find libraries from biopsy")
    parser.add_argument("-t", "--tumour-lib", help="Specify DNA tumour library")
    parser.add_argument("-n", "--normal-lib", help="Specify DNA tumour library")
    parser.add_argument(
        "-o",
        "--output-file",
        help="specify a config filename.",
        default=f"DERIVED_OUTPUT_DIR/{CONFIG_BASENAME}",
    )
    parser.add_argument(
        "-p",
        "--project",
        help="Specify a project instead of bioapps lookup by patient_id.",
    )

    parser.add_argument("-d", "--use-docker", help="Set use docker/slurm tag")
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
        help="Show version and exit.",
    )

    args = parser.parse_args()
    if (args.tumour_lib or args.normal_lib) and args.biopsy:
        raise ValueError("Specify a patient biopsy or libraries, not both.")
    if not any((args.tumour_lib, args.normal_lib, args.biopsy)):
        raise ValueError("Specify a patient biopsy or tumour-lib and normal-lib")
    return args


def main(args=None):
    """Build a GSC config for running Ploidetect.

    Uses GSC standars patient_id and biopsy or libraries.
    Ploidetect versions must just be manually specified.
    """
    args = parse_args() if not args else args

    if args.output_file and exists(args.output_file):
        raise ValueError(f"Output config already exists: '{args.output_file}'")

    config = build_config(**vars(args))

    if not args.output_file:
        YAML().dump(config, sys.stdout)
    else:
        args.output_file = args.output_file.replace(
            "DERIVED_OUTPUT_DIR", config["output_dir"]
        )
        if exists(args.output_file):
            raise ValueError(f"Output config already exists: '{args.output_file}'")
        elif (dirname(args.output_file)) and not exists(dirname(args.output_file)):
            logging.warning(f"Creating output folder: {dirname(args.output_file)}")
            os.makedirs(dirname(args.output_file))
        print(f"Writing config to: {abspath(realpath(args.output_file))}")
        YAML().dump(config, open(args.output_file, "w"))


if __name__ == "__main__":
    main()
