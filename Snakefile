import glob
import os
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import requests

sys.path.insert(0, workflow.basedir)
from constants import VERSION
from constants import WorkflowSetupError

HTTP = HTTPRemoteProvider()
MEM_PER_CPU = 7900

print(f"Ploidetect-pipeline {VERSION}")


## Load config values
configfile: os.path.join(workflow.basedir, "config/defaults.yaml")
configfile: os.path.join(workflow.basedir, "config/samples.yaml")
configfile: os.path.join(workflow.basedir, "config/parameters.yaml")


MEM_PER_CPU = config["mem_per_cpu"] if "mem_per_cpu" in config else MEM_PER_CPU

##### Default checks
## If chromosomes not specified, try to use defaults

if "chromosomes" not in config:
    if config["genome_name"] in config["chromosome_defaults"].keys():
        chromosomes = config["chromosome_defaults"][config["genome_name"]]
    else:
        chromosomes = []
else:
    chromosomes = config["chromosomes"]

if len(chromosomes) == 0:
    raise (
        WorkflowSetupError(
            "No valid chromosomes found. If you're using a non-standard genome (not hg38 or hg19), you must explicitly specify chromosome names in the config"
        )
    )

if "genome" in config:
    if config["genome_name"] in config["genome"]:
        genome_path = config["genome"][config["genome_name"]]
    else:
        genome_path = "failed"
    if not os.path.exists(genome_path):
        genome_path = "failed"
else:
    genome_path == "failed"

if genome_path == "failed":
    hgver = config["genome_name"]
    connec = requests.get(
        f"http://hgdownload.cse.ucsc.edu/goldenpath/{hgver}/database/cytoBand.txt.gz"
    )
    if connec.status_code == 200:
        genome_path = f"resources/{hgver}/genome.fa"
    else:
        raise (
            WorkflowSetupError(
                "Genome fasta file not found in config or on UCSC. Specify a path to your genome fasta"
            )
        )


include: "defaults.smk"


output_dir = config["output_dir"] if "output_dir" in config else "ploidetect_out"

scripts_dir = os.path.join(workflow.basedir, "scripts")
array_positions = (
    config["array_positions"][config["genome_name"]]
    if os.path.exists(config["array_positions"][config["genome_name"]])
    else os.path.join(
        workflow.basedir, config["array_positions"][config["genome_name"]]
    )
)

if "maxd" in config and "qual" in config:
    pass
elif "sequence_type" in config:
    config["maxd"] = config["sequence_type_defaults"][config["sequence_type"]]["maxd"]
    config["qual"] = config["sequence_type_defaults"][config["sequence_type"]]["qual"]
else:
    logger.warning(
        f"No sequence_type or maxd/qual parameters.  Using 'short' defaults."
    )
    config["maxd"] = config["sequence_type_defaults"]["short"]["maxd"]
    config["qual"] = config["sequence_type_defaults"]["short"]["qual"]
if "window_threshold" not in config:
    config["window_threshold"] = config["default_window_threshold"]
if "ploidetect_ver" not in config:
    config["ploidetect_ver"] = config["default_ploidetect_ver"]

with open(array_positions) as f:
    snp_positions = f.readlines()
    snp_positions = [snp.strip().split("\t") for snp in snp_positions]
    snp_chrs = [line[0] for line in snp_positions]
    snp_chrs = list(set(snp_chrs))
    for chromosome in chromosomes:
        if str(chromosome) not in snp_chrs:
            raise (
                WorkflowSetupError(
                    f"Chromosome {chromosome} specified in config was not found in provided snp position data"
                )
            )
    lens = [len(line) for line in snp_positions]
    print(lens[0])
    if any(leng > 2 for leng in lens):
        raise (
            WorkflowSetupError(
                "More than two columns detected in snp position data - file must contain only chromosome and position columns"
            )
        )


cyto_path = config["cyto_path"] if "cyto_path" in config else "auto"
if cyto_path == "auto":
    ## Try to automatically get cytobands based on genome name
    hgver = config["genome_name"]
    connec = requests.get(
        f"http://hgdownload.cse.ucsc.edu/goldenpath/{hgver}/database/cytoBand.txt.gz"
    )
    if connec.status_code == 200:
        cyto_path = f"{workflow.basedir}/resources/{hgver}/cytobands.txt"
    else:
        raise WorkflowSetupError(
            "Cytoband file is set to auto-detect, but could not download cytoband file. Make sure you didn't misspell the genome file, leave the cyto_path blank in the config, or explicitly set a path for it"
        )
cyto_arg = f"-c {cyto_path}" if cyto_path else ""

## Parse sample information
bams_dict = config["bams"]
sample_ids = bams_dict.keys()

## Here comes the for loop
## Initialize output list, which will be fed to rule all
output_list = []
## Loop over the bams and construct lib comparisons
for sample in sample_ids:
    somatics = bams_dict[sample]["somatic"].keys()
    somatic_paths = bams_dict[sample]["somatic"].values()
    normals = bams_dict[sample]["normal"].keys()
    normal_paths = bams_dict[sample]["normal"].values()
    ## Check that the bams are findable
    concat_bams = list(somatic_paths) + list(normal_paths)
    for bam in concat_bams:
        if not os.path.exists(bam):
            raise WorkflowSetupError(
                f"Input file {bam} could not be found. Ensure that the file is spelled correctly, and that you've correctly bound the directory if using singularity"
            )
    combinations = expand("{somatic}_{normal}", somatic=somatics, normal=normals)
    outs = [os.path.join(output_dir, sample, comb, "cna.txt") for comb in combinations]
    output_list.extend(outs)

print(f"Final outputs: {output_list}")


rule all:
    input:
        [
            expand(
                "{output_dir}/{case}/{somatic}_{normal}/cna.txt",
                output_dir=output_dir,
                case=case,
                somatic=config["bams"][case]["somatic"].keys(),
                normal=config["bams"][case]["normal"].keys(),
            )
            for case in config["bams"].keys()
        ],


# Build ploidetect install string
ploidetect_install_cmd = (
    "devtools::install_github('lculibrk/Ploidetect', ref = '{}')".format(
        config["ploidetect_ver"]
    )
)
if "ploidetect_local_clone" in config and config["ploidetect_local_clone"]:
    install_path = config["ploidetect_local_clone"].format(**config)
    ploidetect_install_cmd = f"devtools::install_local('{install_path}', force = TRUE)"
# Bypass untar error - https://github.com/ContinuumIO/anaconda-issues/issues/11509
SET_TAR = "Sys.setenv(TAR = \'/bin/tar\'); "
ploidetect_install_cmd = f'"{SET_TAR + ploidetect_install_cmd}"'


rule ploidetect_install:
    """Install Ploidetect R script into environment"""
    output:
        expand(
            "{install_dir}/conda_configs/ploidetect_installed_{ploidetect_ver}.txt",
            install_dir=workflow.basedir,
            ploidetect_ver=config["ploidetect_ver"],
        ),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    params:
        install_cmd=ploidetect_install_cmd,
    log:
        expand("{output_dir}/logs/ploidetect_install.log", output_dir=output_dir),
    shell:
        "date >> {log}; echo {params} >> {log}; "
        "export LC_ALL=en_US.UTF-8; "
        " Rscript -e {params.install_cmd} 2>> {log}"
        " && echo {params} > {output} && date >> {output}"


rule download_cytobands:
    """Downloads cytoband data for plotting & (todo) hgver-specific centromere filtering"""
    input:
        HTTP.remote(
            expand(
                "http://hgdownload.cse.ucsc.edu/goldenpath/{hgver}/database/cytoBand.txt.gz",
                hgver=config["genome_name"],
            )
        ),
    output:
        expand(
            "{install_dir}/resources/{hgver}/cytobands.txt",
            install_dir=workflow.basedir,
            hgver=config["genome_name"],
        ),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    shell:
        "gunzip -c {input} | sed 's/chr//g' > {output}"


rule download_genome:
    input:
        HTTP.remote(
            expand(
                "https://hgdownload.soe.ucsc.edu/goldenPath/{hgver}/bigZips/{hgver}.fa.gz",
                hgver=config["genome_name"],
            )
        ),
    output:
        expand("resources/{hgver}/genome.fa", hgver=config["genome_name"]),
        expand("resources/{hgver}/genome.fa.fai", hgver=config["genome_name"]),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "gunzip -c {input} > {output[0]}; samtools faidx {output[0]}"


rule germline_cov:
    """Compute per-base depth in germline bam, convert to .bed format and pile up into equal-coverage bins"""
    input:
        bam=lambda w: config["bams"][w.case]["normal"][w.normal],
    output:
        temp("{output_dir}/scratch/{case}/{normal}/normal/{chr}.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
        threshold=config["window_threshold"],
        qual=config["qual"],
        maxd=config["maxd"],
    log:
        "{output_dir}/logs/germline_cov.{case}.{normal}.{chr}.log",
    shell:
        "samtools depth -r{wildcards.chr} -Q{params.qual} -m {params.maxd} {input.bam} 2>> {log}"
        " | awk -v FS='\\t' -v OFS='\\t' 'NR > 1{{print $1, $2, $2+1, $3}}'"
        " | python3 {params.scripts_dir}/make_windows.py - {params.threshold} 2>> {log}"
        " | bedtools sort -i stdin > {output}  2>> {log}"


rule merge_germline:
    """Merge multi-chromosome output from germline_cov into single file"""
    input:
        expand(
            "{{output_dir}}/scratch/{{case}}/{{normal}}/normal/{chr}.bed",
            chr=chromosomes,
        ),
    output:
        temp("{output_dir}/scratch/{case}/{normal}/germline.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/merge_germline.{case}.{normal}.log",
    shell:
        "ls -l {input} >> {log};"
        "cat {input} | bedtools sort -i stdin > {output}"
        " 2>> {log}"


rule makewindowfile:
    """Remove germline depth column from file to obtain bins"""
    input:
        rules.merge_germline.output,
    output:
        temp("{output_dir}/scratch/{case}/{normal}/windows.txt"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/makewindowfile.{case}.{normal}.log",
    shell:
        "cut -f 1,2,3 < {input} | bedtools sort -i stdin > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule splitwindowfile:
    """Split bins into each chromosome for parallel computing of depth in somatic"""
    input:
        rules.makewindowfile.output,
    output:
        temp("{output_dir}/scratch/{case}/{normal}/windows/{chr}.txt"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/splitwindowfile.{case}.{normal}.{chr}.log",
    shell:
        "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule genomecovsomatic:
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        window=rules.splitwindowfile.output,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/tumour/{chr}.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    params:
        qual=config["qual"],
        maxd=config["maxd"],
    log:
        "{output_dir}/logs/genomecovsomatic.{case}.{somatic}_{normal}.{chr}.log",
    shell:
        "samtools depth -Q {params.qual} -m {params.maxd} -r {wildcards.chr} -a {input.sombam} "
        " | awk -v FS='\\t' -v OFS='\\t' '{{print $1, $2, $2 + 1, $3}}'"
        " | sort -k1,1 -k2,2n "
        " | bedtools map -b stdin -a {input.window} -c 4 -o mean > {output}"


rule genomecovgermline:
    input:
        normbam=lambda w: config["bams"][w.case]["normal"][w.normal],
        window=rules.splitwindowfile.output,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/normal/{chr}.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    params:
        qual=config["qual"],
        maxd=config["maxd"],
    log:
        "{output_dir}/logs/genomecovgermline.{case}.{somatic}_{normal}.{chr}.log",
    shell:
        "samtools depth -Q {params.qual} -m {params.maxd} -r {wildcards.chr} -a {input.normbam}"
        " | awk -v FS='\\t' -v OFS='\\t' '{{print $1, $2, $2 + 1, $3}}'"
        " | sort -k1,1 -k2,2n"
        " | bedtools map -b stdin -a {input.window} -c 4 -o mean > {output}"


rule merge_split_tumour:
    input:
        expand(
            "{{output_dir}}/scratch/{{case}}/{{somatic}}_{{normal}}/tumour/{chr}.bed",
            chr=chromosomes,
        ),
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/tumour.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    log:
        "{output_dir}/logs/mergesomatic.{case}.{somatic}_{normal}.log",
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule merge_split_normal:
    input:
        expand(
            "{{output_dir}}/scratch/{{case}}/{{somatic}}_{{normal}}/normal/{chr}.bed",
            chr=chromosomes,
        ),
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/normal.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    log:
        "{output_dir}/logs/mergenormal.{case}.{somatic}_{normal}.log",
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule compute_loh:
    """Variant call the germline, filter for heterozygous snps and count alleles in somatic"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        normbam=lambda w: config["bams"][w.case]["normal"][w.normal],
        genome=genome_path,
    output:
        folder=directory("{output_dir}/scratch/{case}/{somatic}_{normal}/loh_tmp"),
        loh=temp("{output_dir}/scratch/{case}/{somatic}_{normal}/loh_tmp/loh_raw.txt"),
    params:
        array_positions={array_positions},
        scripts_dir=scripts_dir,
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/compute_loh.{case}.{somatic}_{normal}.log",
    shell:
        "bash {params.scripts_dir}/get_allele_freqs.bash {input.normbam} {input.sombam}"
        " {input.genome} {params.array_positions}"
        " {output.folder}"
        " &> {log}"


rule process_loh:
    """Convert allele counts to beta-allele frequencies and merge for each bin"""
    input:
        loh=rules.compute_loh.output.loh,
        window=rules.makewindowfile.output,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/loh.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
    log:
        "{output_dir}/logs/process_loh.{case}.{somatic}_{normal}.log",
    shell:
        "awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}' {input.loh}"
        " | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}'"
        " | bedtools sort -i stdin"
        " | Rscript {params.scripts_dir}/merge_loh.R -l STDIN -w {input.window} -o {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule getgc:
    """Get GC content for each bin"""
    input:
        window=rules.makewindowfile.output,
        genome=genome_path,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/gc.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/getgc.{case}.{somatic}_{normal}.log",
    shell:
        "bedtools nuc -fi {input.genome} -bed {input.window} | cut -f1,2,3,5 | tail -n +2 > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule mergedbed:
    """Merge all the data into a single file.  Remove chr prefix from chromsomes."""
    input:
        gc=rules.getgc.output,
        tumour="{output_dir}/scratch/{case}/{somatic}_{normal}/tumour.bed",
        normal="{output_dir}/scratch/{case}/{somatic}_{normal}/normal.bed",
        loh=rules.process_loh.output,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/merged.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/mergedbed.{case}.{somatic}_{normal}.log",
    shell:
        "paste {input.tumour} <(cut -f4 {input.normal}) <(cut -f4 {input.loh}) <(cut -f4 {input.gc})"
        "| sed 's/chr//g' > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        merged="{output_dir}/scratch/{case}/{somatic}_{normal}/merged.bed",
        install=rules.ploidetect_install.output
        if workflow.use_conda and not workflow.use_singularity
        else __file__,
        cytos=cyto_path,
    output:
        "{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
    log:
        "{output_dir}/logs/preseg.{case}.{somatic}_{normal}.log",
    shell:
        "Rscript {params.scripts_dir}/prep_ploidetect2.R -i {input.merged} -c {input.cytos} -o {output}"


rule ploidetect:
    """Runs Ploidetect"""
    input:
        preseg=rules.preseg.output,
        install=(
            rules.ploidetect_install.output
            if not workflow.use_singularity
            and "install_ploidetect" in config.keys()
            and config["install_ploidetect"]
            else __file__
        ),
    output:
        plots="{output_dir}/{case}/{somatic}_{normal}/plots.pdf",
        models="{output_dir}/{case}/{somatic}_{normal}/models.txt",
        meta="{output_dir}/{case}/{somatic}_{normal}/meta.RDS",
    conda:
        "conda_configs/r.yaml"
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
        cyto_arg=cyto_arg,
    log:
        "{output_dir}/logs/ploidetect.{case}.{somatic}_{normal}.log",
    shell:
        "Rscript {params.scripts_dir}/run_ploidetect2.R "
        " -i {input.preseg} "
        " -o {output.models} -p {output.plots} -r {output.meta} {params.cyto_arg}"
        " &> {log}"


rule force_tcp:
    """Forces a new model of purity/ploidy for CNV calling if specified in the config"""
    input:
        model="{output_dir}/{case}/{somatic}_{normal}/models.txt",
    output:
        model="{output_dir}/{case}/{somatic}_{normal}/models_cnv.txt",
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=1 * MEM_PER_CPU,
    params:
        tp=(
            lambda w: config["models"][w.case][w.somatic]["tp"]
            if "models" in config
            and w.case in config["models"]
            and w.somatic in config["models"][w.case]
            else "NA"
        ),
        ploidy=(
            lambda w: config["models"][w.case][w.somatic]["ploidy"]
            if "models" in config
            and w.case in config["models"]
            and w.somatic in config["models"][w.case]
            else "NA"
        ),
    ## CNV caller's log to record if a non-automated tc/ploidy was used
    log:
        "{output_dir}/logs/ploidetect_copynumber.{case}.{somatic}_{normal}.log",
    shell:
        "Rscript -e '\n "
        "require(data.table) \n"
        'd = suppressWarnings(fread("{input}")) \n '
        "   if(is.na({params.tp})){{ \n "
        '       message("CNV calling using automatically detected purity/ploidy values") \n'
        '       fwrite(d, "{output}", sep = "\\t") \n '
        "       quit(status = 0) \n "
        "   }} \n "
        '   message(paste0("Manually provided purity of ", {params.tp}, " and ploidy of ", {params.ploidy}, " specified, using those.")) \n'
        "   d$tp[1] = {params.tp} \n "
        "   d$ploidy[1] = {params.ploidy} \n "
        '   fwrite(d, "{output}", sep = "\\t")\' 2> {log}'


rule ploidetect_copynumber:
    """Performs CNV calling using the tumor purity and ploidy estimated by Ploidetect"""
    input:
        cytos=cyto_path,
        models="{output_dir}/{case}/{somatic}_{normal}/models.txt",
        rds="{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
        plots="{output_dir}/{case}/{somatic}_{normal}/plots.pdf",
    output:
        cna="{output_dir}/{case}/{somatic}_{normal}/cna.txt",
        cna_plots="{output_dir}/{case}/{somatic}_{normal}/cna_plots.pdf",
        cna_cond="{output_dir}/{case}/{somatic}_{normal}/cna_condensed.txt",
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    params:
        scripts_dir=scripts_dir,
        cyto_arg=cyto_arg,
    log:
        "{output_dir}/logs/ploidetect_copynumber.{case}.{somatic}_{normal}.log",
    shell:
        "Rscript {params.scripts_dir}/ploidetect_copynumber.R"
        " -i {input.rds} -m {input.models}"
        " -p {output.cna_plots} -o {output.cna} {params.cyto_arg}"
        " &> {log}"
