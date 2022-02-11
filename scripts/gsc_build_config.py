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

__version__ = "0.1.0"

logger = logging.getLogger(__name__)
API = BioappsApi()
CONFIG_BASENAME = "Ploidetect-pipeline.yaml"

REFS = ("hg19", "hg38")
DEFAULTS_FOLDER = abspath(join(dirname(abspath(__file__)), "..", "resources", "config"))
GENOME_REF_YAMLS = [join(DEFAULTS_FOLDER, f"genome_ref.{ref}.yaml") for ref in REFS]
DEFAULT_RUN_YAML = join(DEFAULTS_FOLDER, "default_run_params.yaml")


def get_biopsy_dna_tumour_normal(id, biopsy="biop1"):
    """Return DNA (tumour, normal) for the patient biopsy.

    Example:
        >>> get_biopsy_dna_tumour_normal("POG965", biopsy="biop2")
        ('P02866', 'P02590')
    """
    libs = API.get_lib_names(id)
    normals = [
        lib for lib in libs if API.is_normal_library(lib) and API.is_dna_library(lib)
    ]
    tumour_libs = [
        lib for lib in libs if API.get_biopsy(lib) == biopsy and API.is_dna_library(lib)
    ]
    assert len(normals) == 1, f"Multiple normal libs: {normals}"
    assert len(tumour_libs) == 1, f"Multiple {biopsy} tumour_libs: {tumour_libs}"
    return (tumour_libs[0], normals[0])


def get_project_info(id):
    """Return a single project_info dict for the patient."""
    proj_info = API.get_biofx_project_info(id, tumour_char_only=True)
    if len(proj_info) > 1:
        logger.error(
            f"Assuming project: {proj_info[0]['name']}"
            f" - ignoring {[pi['name'] for pi in proj_info[1:]]}"
        )
    return proj_info[0]


def get_somatic_folder(id, tumour_lib, normal_lib, project=None, genome_name=None):
    """Return the 'somatic' paired analysis folder.

    Example:
        >>> tumour,normal = get_biopsy_dna_tumour_normal("POG965", biopsy="biop2")
        >>> get_somatic_folder("POG965", tumour, normal)
        '/projects/POG/POG_data/POG965/wgs/biop2_t_P02866_blood1_n_P02590'
    """
    patient_proj_info = get_project_info(id)
    patient_project = patient_proj_info["name"]
    patient_project_path = patient_proj_info["path"]
    project = project if project else patient_project

    if patient_project != project:
        logger.error(f"Project conflict - {project} vs {id} defined {patient_project}")

    project_path = API.get_biofx_project_path(project, genome_name) or patient_project_path
    if patient_project_path != project_path:
        logger.error(
            f"Conflicting project paths - using ({project}, {genome_name}): {project_path} instead of {id}: {patient_project_path}"
        )

    glob_str = f"{project_path}/{id}/*/*_t_{tumour_lib}*_n_{normal_lib}"
    somatic_folders = glob.glob(glob_str)
    if not somatic_folders:
        raise ValueError(f"Failed to find somatic folder.  Checked '{glob_str}'")
    elif len(somatic_folders) > 1:
        raise ValueError(
            f"Failed to find a unique somatic folder.  Found {somatic_folders}"
        )
    return somatic_folders[0]


def get_ploidetect_temp_folder(
    id, tumour_lib, normal_lib, pipeline_ver, ploidetect_ver, project=None, genome_name=None, job_id=None
):
    """GSC trans_scratch temp folder.

    Example:
        >>> get_ploidetect_temp_folder("COLO829-TestA", "A36971", "A36973", "0.0.1", "v1.0.0", "POG")
        '/projects/trans_scratch/validations/Ploidetect/POG/COLO829-TestA/Ploidetect-pipeline-0.0.1/Ploidetect-v1.0.0/A36971_A36973'
    """
    TEMP_ROOT = "/projects/trans_scratch/validations/Ploidetect"
    if not project:
        project = get_project_info(id)["name"]
    scratch_path = join(
        TEMP_ROOT,
        project,
        id,
        "Ploidetect-pipeline-" + pipeline_ver,
        "Ploidetect-" + ploidetect_ver,
        f"{tumour_lib}_{normal_lib}",
    )
    if genome_name:
        scratch_path = join(scratch_path, genome_name)
    if job_id:
        scratch_path = join(scratch_path, job_id)
    return scratch_path


def get_gsc_output_folder(
    id,
    tumour_lib,
    normal_lib,
    pipeline_ver,
    ploidetect_ver,
    project=None,
    genome_name=None,
):
    """GSC final project output folder.

    Example:
        >>> get_gsc_output_folder("POG965", "P02866", "P02590", "0.0.1", "v1.0.0")
        '/projects/POG/POG_data/POG965/wgs/biop2_t_P02866_blood1_n_P02590/Ploidetect/Ploidetect-pipeline-0.0.1/Ploidetect-v1.0.0'
    """
    somatic_pair_folder = get_somatic_folder(
        id, tumour_lib, normal_lib, project=project, genome_name=genome_name
    )
    retval = join(
        somatic_pair_folder,
        "Ploidetect",
        "Ploidetect-pipeline-" + pipeline_ver,
        "Ploidetect-" + ploidetect_ver,
    )
    if genome_name:
        retval = join(retval, genome_name)
    return retval


def get_bam(library, genome_name=None):
    """Find a bam for library.  Restricted genome_refernce if given.

    Returns (library, genome_reference)
    """
    bams_all = API.get_bam_info_from_library(library)
    if genome_name:
        bams_all = [
            bd for bd in bams_all if bd["genome_reference"].startswith(genome_name)
        ]

    bams = []
    for bam in bams_all:
        bam_fns = glob.glob(join(bam["data_path"], "*.bam"))
        if not bam_fns:  # GERO-114 - COLO829 - no bam in path
            logger.error(f"No bam found in path: {bam['data_path']}")
        else:
            bams.append(bam)

    assert bams, f"Bam finding error for: '{library}' (ref: {genome_name})"

    bams.sort(key=lambda x: bool(x["tumour_char"]), reverse=True)
    for bam in bams[1:]:
        warn = f"Ignoring tumour bam {bam['data_path']}"
        logger.warning(warn)

    bam_fns = glob.glob(join(bams[0]["data_path"], "*.bam"))
    assert (
        len(bam_fns) == 1
    ), f"Bam finding error for: '{library}' (ref: {genome_name})"
    return (bam_fns[0], bams[0]["genome_reference"])


def genome_reference2genome_name(gsc_genome_reference):
    if gsc_genome_reference.lower() in ("hg19", "hg19a", "grch37"):
        return "hg19"
    if gsc_genome_reference.lower() in ("hg38", "hg38_no_alt", "grch38"):
        return "hg38"
    return gsc_genome_reference


def build_config(
    id,
    biopsy=None,
    tumour_lib=None,
    normal_lib=None,
    genome_name=None,
    pipeline_ver="undefined",
    ploidetect_ver="undefined",
    install_ploidetect=False,
    project=None,
    output_dir="",
    ref_yamls=GENOME_REF_YAMLS,
    run_params_yaml=DEFAULT_RUN_YAML,
    **kwargs,
):
    """Build a GSC config for running Ploidetect.

    Checks GSC bioapps api for relevant details of libraries, bams and output paths.

    Example:
        >>> cfg = build_config("POG965", biopsy="biop2")
        >>> YAML().dump(cfg, sys.stdout)  # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
        # Created by: ...
        bams:
          POG965:
            somatic:
              P02866: /projects/analysis/analysis30/P02866/merge32550_bwa-mem-0.7.6a-sb/150bp/hg19a/P02866_3_lanes_dupsFlagged.bam
            normal:
              P02590: /projects/analysis/analysis30/P02590/HCW32CCXY_8/P02590/150nt/hg19a/bwa-mem-0.7.6a-sb/P02590_1_lane_dupsFlagged.bam
        genome_name: hg19
        output_dir: /projects/POG/POG_data/POG965/wgs/biop2_t_P02866_blood1_n_P02590/Ploidetect/Ploidetect-pipeline-undefined/Ploidetect-undefined
        temp_dir: /projects/trans_scratch/validations/Ploidetect/POG/POG965/Ploidetect-pipeline-undefined/Ploidetect-undefined/P02866_P02590...
        # ploidetect_ver should be a branch or tag.  Overriden by ploidetect_local_clone.
        ploidetect_ver: undefined
        # Leave ploidetect_local_clone blank or 'None' to download from github
        ploidetect_local_clone: /gsc/pipelines/Ploidetect/undefined
        install_ploidetect: 0
        window_threshold: 100000
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
        - Y
    """
    TAB = "  "
    yaml = YAML()

    # Keep input parameters
    prog_str = f'Created by: {realpath(abspath(__file__))} at {datetime.now().strftime("%Y%m%d %H:%M:%S")}'
    yaml_lines = [f"# {prog_str}"]
    yaml_lines.append(f"id: {id}")
    if biopsy:
        yaml_lines.append(f"biopsy: {biopsy}")

    if not biopsy and not (tumour_lib and normal_lib):
        raise ValueError(
            "Either a biopsy or tumour_lib and normal_lib must be supplied."
        )
    elif biopsy:
        tumour_lib, normal_lib = get_biopsy_dna_tumour_normal(id, biopsy)

    yaml_lines.append(f"tumour_lib: {tumour_lib}")
    yaml_lines.append(f"normal_lib: {normal_lib}")

    # Find bams
    yaml_lines.append("bams:")
    yaml_lines.append(f"{TAB}{id}:")
    tumour_bam_fn, genome_name = get_bam(tumour_lib, genome_name)
    normal_bam_fn, normal_genome_name = get_bam(normal_lib, genome_name)

    yaml_lines.append(f"{TAB}{TAB}somatic:")
    yaml_lines.append(f"{TAB}{TAB}{TAB}{tumour_lib}: {tumour_bam_fn}")
    yaml_lines.append(f"{TAB}{TAB}normal:")
    yaml_lines.append(f"{TAB}{TAB}{TAB}{normal_lib}: {normal_bam_fn}")

    assert genome_name == normal_genome_name
    yaml_lines.append(f"genome_name: {genome_reference2genome_name(genome_name)}")

    if not output_dir:
        output_dir = get_gsc_output_folder(id,
                                           tumour_lib,
                                           normal_lib,
                                           pipeline_ver,
                                           ploidetect_ver,
                                           project=project,
                                           genome_name=genome_name)
    yaml_lines.append(f"output_dir: {output_dir}")

    # trans_scratch temp directory
    if "temp_dir" in kwargs:
        temp_dir = kwargs["temp_dir"]
    else:
        temp_dir = get_ploidetect_temp_folder(
            id,
            tumour_lib,
            normal_lib,
            pipeline_ver,
            ploidetect_ver,
            project=project,
        )
        temp_dir = join(temp_dir, datetime.now().strftime("%Y%m%d_%H%M%S"))
        logger.warning(f"Derived GSC scratch folder: {temp_dir}")
    yaml_lines.append(f"temp_dir: {temp_dir}")

    # Ploidetect installation and versions.
    yaml_lines.append(
        "# ploidetect_ver should be a branch or tag.  Overriden by ploidetect_local_clone."
    )
    yaml_lines.append(f"ploidetect_ver: {ploidetect_ver}")
    yaml_lines.append(
        "# Leave ploidetect_local_clone blank or 'None' to download from github"
    )
    yaml_lines.append(
        f"ploidetect_local_clone: /gsc/pipelines/Ploidetect/{ploidetect_ver}"
    )
    yaml_lines.append(f"install_ploidetect: {1 if bool(install_ploidetect) else 0}")

    config = yaml.load("\n".join(yaml_lines))
    # Store all additional values given
    for k, v in kwargs.items():
        if k not in config.keys():
            config[k] = v

    logger.info(f"Loading genome reference data from f{run_params_yaml}")
    run_params = yaml.load(open(run_params_yaml))
    for k, v in run_params.items():
        if k not in config.keys():
            logger.warning(f"Adding {k} from {run_params_yaml}")
            config[k] = v

    # Genomic Reference data
    # each reference should have the same keys:
    #   ['annotation', 'array_positions', 'genome', 'ref_chromosomes']
    # Values for reference should be under reference names to avoid conflicts. eg:
    #   ref['annotation']['hg38'] = <filepath_to_gtf>
    for ref_yaml_fn in ref_yamls:
        logger.info(f"Loading genome reference data from {ref_yaml_fn}")
        ref = yaml.load(open(ref_yaml_fn))
        for k, v in ref.items():
            if k not in config:
                config[k] = v
            else:
                config[k].update(v)
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
        "--id",
        required=True,
        help="GSC bioapps patient_id.  eg POG965",
    )
    parser.add_argument("-b", "--biopsy", help="Find libraries from biopsy")
    parser.add_argument("-t", "--tumour-lib", help="Specify DNA tumour library")
    parser.add_argument("-n", "--normal-lib", help="Specify DNA tumour library")
    parser.add_argument(
        "-r", "--genome-reference", help="Specify hg19/hg38 instead of using tc flag"
    )
    parser.add_argument(
        "-o",
        "--gsc-config-filename",
        help="specify a config filename.",
        default=f"DERIVED_OUTPUT_DIR/{CONFIG_BASENAME}",
    )
    parser.add_argument("--output_dir", help="Output directory override.")
    parser.add_argument(
        "-p",
        "--project",
        help="Specify a project instead of bioapps lookup by patient_id.",
    )

    parser.add_argument("--install_ploidetect", help="Install into local R.")
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
    logger.setLevel(logging.INFO)

    args = parse_args() if not args else args

    if args.gsc_config_filename and exists(args.gsc_config_filename):
        raise ValueError(f"Output config already exists: '{args.gsc_config_filename}'")

    config = build_config(**vars(args))

    if not args.gsc_config_filename:
        YAML().dump(config, sys.stdout)
    else:
        args.gsc_config_filename = args.gsc_config_filename.replace(
            "DERIVED_OUTPUT_DIR", config["output_dir"]
        )
        if exists(args.gsc_config_filename):
            raise ValueError(f"Output config already exists: '{args.gsc_config_filename}'")
        elif (dirname(args.gsc_config_filename)) and not exists(dirname(args.gsc_config_filename)):
            logger.warning(f"Creating output folder: {dirname(args.gsc_config_filename)}")
            os.makedirs(dirname(args.gsc_config_filename))
        print(f"Writing config to: {abspath(realpath(args.gsc_config_filename))}")
        YAML().dump(config, open(args.gsc_config_filename, "w"))


if __name__ == "__main__":
    main()
