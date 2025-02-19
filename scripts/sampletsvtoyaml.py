#!/usr/bin/env python
"""sampletsvtoyaml.py: Converts tab-separated sample information to configuration file for Ploidetect"""

import argparse
import logging
import os

import yaml


class inputFileError(ValueError):
    """raise this when there's a mistake in an input file"""


parser = argparse.ArgumentParser(
    description="Converts a tab-separated file of sample information to a samples.yaml file for the Ploidetect pipeline"
)

parser.add_argument("-t", "--tsv_file", help="Path to input tsv file")
parser.add_argument(
    "-o",
    "--output_path",
    help="Path to config directory. Defaults to ./config/",
    default="./config/",
)
parser.add_argument(
    "-c",
    "--concat",
    help="Should the output be concatenated onto an existing samples.yaml? Accepts True or False",
)

args = parser.parse_args()

if not os.path.exists(args.output_path):
    os.mkdir(args.output_path)

with open(args.tsv_file) as f:
    file_rows = [line.strip().split("\t") for line in f.readlines()]

lineno = 0
out_dict = {"bams": {}}
for row in file_rows:
    samplename = row[0]
    identity = row[1]
    lib = row[2]
    path = row[3]
    lineno = lineno + 1
    if len(row) != 4:
        raise inputFileError(
            f"Input samples file does not have four tab-separated columns! Error detected in line {lineno}"
        )
    if not os.path.exists(path):
        logging.warning(
            f"Input samples file points to a file which does not exist. Problematic path: {path} on line {lineno}. If this is on purpose, ignore this warning"
        )
    if not bool(out_dict["bams"]):
        out_dict["bams"][samplename] = {}
    out_dict["bams"][samplename][identity] = {lib: path}

if args.concat == "True":
    out_dict = out_dict["bams"]
    with open(os.path.join(args.output_path, "samples.yaml"), "r") as f:
        cur_dict = yaml.load(f, Loader=yaml.FullLoader)
        n_dict = {**cur_dict["bams"], **out_dict}
        n_dict = {"bams": n_dict}
    with open(os.path.join(args.output_path, "samples.yaml"), "w") as f:
        yaml.dump(n_dict, f, default_flow_style=False)
else:
    with open(os.path.join(args.output_path, "samples.yaml"), "w") as f:
        yaml.dump(out_dict, f, default_flow_style=False)
