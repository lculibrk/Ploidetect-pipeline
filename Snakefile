import glob
import os
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import requests
HTTP = HTTPRemoteProvider()

sys.path.insert(0, workflow.basedir)
from constants import VERSION
from constants import WorkflowSetupError

print(f"Ploidetect-pipeline {VERSION}")


## Load config values
configfile: os.path.join(workflow.basedir, "defaults.yaml")
configfile: os.path.join(workflow.basedir, "CONFIG.txt")

MEM_PER_CPU = config["mem_per_cpu"]

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
    raise(WorkflowSetupError("No valid chromosomes found. If you're using a non-standard genome (not hg38 or hg19), you must explicitly specify chromosome names in the config"))

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
    connec = requests.get(f"http://hgdownload.cse.ucsc.edu/goldenpath/{hgver}/database/cytoBand.txt.gz")
    if connec.status_code == 200:
        genome_path = f"resources/{hgver}/genome.fa"
    else:
        raise(WorkflowSetupError("Genome fasta file not found in config or on UCSC. Specify a path to your genome fasta"))

output_dir = config["output_dir"]
if "temp_dir" not in config or not config["temp_dir"]:
    config["temp_dir"] = f"{output_dir}/temp"
temp_dir = config["temp_dir"]
if temp_dir[-1] != "/":
    temp_dir += "/"  # Prevents strange case wild-card error.

scripts_dir = os.path.join(workflow.basedir, "scripts")
array_positions = (
    config["array_positions"][config["genome_name"]]
    if os.path.exists(config["array_positions"][config["genome_name"]])
    else os.path.join(
        workflow.basedir, config["array_positions"][config["genome_name"]]
    )
)

cyto_path = config["cyto_path"] if "cyto_path" in config else ""
if cyto_path == "auto":
    ## Try to automatically get cytobands based on genome name
    hgver = config["genome_name"]
    connec = requests.get(f"http://hgdownload.cse.ucsc.edu/goldenpath/{hgver}/database/cytoBand.txt.gz")
    if connec.status_code == 200:
        cyto_path = f"resources/{hgver}/cytobands.txt"
        cyto_arg = f"-c {cyto_path}"
    else:
        raise WorkflowSetupError("Cytoband file is set to auto-detect, but could not download cytoband file. Make sure you didn't misspell the genome file, leave the cyto_path blank in the config, or explicitly set a path for it")


## Parse sample information
bams_dict = config["bams"]
sample_ids = bams_dict.keys()

## Here comes the for loop
## Initialize output list, which will be fed to rule all
output_list = []
## Loop over the bams and construct lib comparisons
for sample in sample_ids:
    somatics = bams_dict[sample]["somatic"].keys()
    normals = bams_dict[sample]["normal"].keys()
    combinations = expand("{somatic}_{normal}", somatic=somatics, normal=normals)
    outs = [os.path.join(output_dir, sample, comb, "cna.txt") for comb in combinations]
    output_list.extend(outs)


print(f"Final outputs: {output_list}")

rule all:
    input:
        output_list,


def devtools_install():
    if config["ploidetect_local_clone"] and config["ploidetect_local_clone"] != "None":
        install_path = config["ploidetect_local_clone"].format(**config)
        devtools_cmd = "\"devtools::install_local('" + install_path + "')\""
    else:
        devtools_cmd = "\"devtools::install_github('lculibrk/Ploidetect', "
        devtools_cmd += "ref = '" + config["ploidetect_ver"] + "')\""
    return devtools_cmd


rule ploidetect_install:
    """Install Ploidetect R script into environment"""
    output:
        expand(
            "{install_dir}/conda_configs/ploidetect_installed.txt",
            install_dir=workflow.basedir,
        ),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    params:
        devtools_install(),
    version:
        config["ploidetect_ver"]
    shell:
        "export LC_ALL=en_US.UTF-8; "
        " Rscript -e {params} "
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
        expand("resources/{hgver}/cytobands.txt", hgver=config["genome_name"]),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    shell:
        "gunzip -c {input} > {output}"

rule download_genome:
    input:
        HTTP.remote(
            expand(
                "https://hgdownload.soe.ucsc.edu/goldenPath/{hgver}/bigZips/{hgver}.fa.gz",
                hgver=config["genome_name"],
            )
        )
    output:
        expand("resources/{hgver}/genome.fa", hgver = config["genome_name"]),
        expand("resources/{hgver}/genome.fa.fai", hgver = config["genome_name"]),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "gunzip -c {input} > {output[0]}; samtools faidx {output[0]}"

rule germline_cov:
    """Compute per-base depth in germline bam, convert to .bed format and pile up into equal-coverage bins"""
    input:
        bam=lambda w: config["bams"][w.case]["normal"][w.lib],
    output:
        temp("{temp_dir}/{case}/{lib}/normal/{chr}.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "samtools depth -r{wildcards.chr} -Q 21 {input.bam}"
        " | awk -v FS='\t' -v OFS='\t' 'NR > 1{{print $1, $2, $2+1, $3}}'"
        " | python3 {scripts_dir}/make_windows.py - 100000"
        " | bedtools sort -i stdin > {output}"


rule merge_germline:
    """Merge multi-chromosome output from germline_cov into single file"""
    input:
        expand("{{temp_dir}}/{{case}}/{{normal}}/normal/{chr}.bed", chr=chromosomes),
    output:
        temp("{temp_dir}/{case}/{normal}/germline.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"


rule makewindowfile:
    """Remove germline depth column from file to obtain bins"""
    input:
        rules.merge_germline.output,
    output:
        temp("{temp_dir}/{case}/{normal}/windows.txt"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "cut -f 1,2,3 < {input} | bedtools sort -i stdin > {output}"


rule splitwindowfile:
    """Split bins into each chromosome for parallel computing of depth in somatic"""
    input:
        rules.makewindowfile.output,
    output:
        temp("{temp_dir}/{case}/{normal}/windows/{chr}.txt"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"


rule genomecovsomatic:
    """Compute depth of tumor and normal reads in previously created bins"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        nombam=lambda w: config["bams"][w.case]["normal"][w.normal],
        window=rules.splitwindowfile.output,
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/tumour/{chr}.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "bedtools multicov -bams {input.sombam} {input.nombam} -q 20 -bed {input.window}  > {output}"


rule mergesomatic:
    """Merge output of genomecovsomatic to a singular file"""
    input:
        expand(
            "{{temp_dir}}/{{case}}/{{somatic}}_{{normal}}/tumour/{chr}.bed",
            chr=chromosomes,
        ),
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/tumour.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"


rule compute_loh:
    """Variant call the germline, filter for heterozygous snps and count alleles in somatic"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        normbam=lambda w: config["bams"][w.case]["normal"][w.normal],
        genome = genome_path
    output:
        temp(directory("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp/")),
        temp("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp/loh_raw.txt"),
    params:
        array_positions={array_positions},
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "bash {scripts_dir}/get_allele_freqs.bash {input.normbam} {input.sombam}"
        " {input.genome} {params.array_positions}"
        " {output[0]}"


rule process_loh:
    """Convert allele counts to beta-allele frequencies and merge for each bin"""
    input:
        loh=rules.compute_loh.output,
        window=rules.makewindowfile.output,
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/loh.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}' {input.loh[1]}"
        " | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}'"
        " | bedtools sort -i stdin"
        " | Rscript {scripts_dir}/merge_loh.R -l STDIN -w {input.window} -o {output}"


rule getgc:
    """Get GC content for each bin"""
    input:
        window=rules.makewindowfile.output,
        genome = genome_path,
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/gc.bed"),
    params:
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "bedtools nuc -fi {input.genome} -bed {input.window} | cut -f1,2,3,5 | tail -n +2 > {output}"


rule mergedbed:
    """Merge all the data into a single file"""
    input:
        gc=rules.getgc.output,
        tumour=rules.mergesomatic.output,
        normal=rules.merge_germline.output,
        loh=rules.process_loh.output,
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/merged.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc}) > {output}"


rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        temp_dir + "{case}/{somatic}_{normal}/merged.bed",
        rules.ploidetect_install.output if not workflow.use_singularity and "install_ploidetect" in config.keys() and config[
            "install_ploidetect"
        ] else __file__,
        cytos=cyto_path,
    output:
        "{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect:devel"
    shell:
        "Rscript {scripts_dir}/prep_ploidetect2.R -i {input[0]} -c {input.cytos} -o {output}"


rule ploidetect:
    """Runs Ploidetect"""
    input:
        rules.preseg.output,
        rules.ploidetect_install.output if not workflow.use_singularity and "install_ploidetect" in config.keys() and config[
            "install_ploidetect"
        ] else __file__,
        cytos=cyto_path,
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
        "docker://lculibrk/ploidetect:devel"
    params:
        cyto_arg = cyto_arg
    shell:
        "Rscript {scripts_dir}/run_ploidetect2.R "
        " -i {input[0]} "
        " -o {output.models} -p {output.plots} -r {output.meta} {params.cyto_arg}"


rule ploidetect_copynumber:
    """Performs CNV calling using the tumor purity and ploidy estimated by Ploidetect"""
    input:
        cytos=cyto_path,
        models="{output_dir}/{case}/{somatic}_{normal}/models.txt",
        rds="{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
        plots="{output_dir}/{case}/{somatic}_{normal}/plots.pdf",
    output:
        "{output_dir}/{case}/{somatic}_{normal}/cna.txt",
        "{output_dir}/{case}/{somatic}_{normal}/cna_plots.pdf",
        "{output_dir}/{case}/{somatic}_{normal}/cna_condensed.txt",
    conda:
        "conda_configs/r.yaml"
    log:
        "{output_dir}/{case}/{somatic}_{normal}/cna_log.txt",
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    container:
        "docker://lculibrk/ploidetect:devel"
    params:
        cyto_arg = cyto_arg
    shell:
        "Rscript {scripts_dir}/ploidetect_copynumber.R -i {input.rds} -m {input.models} -p {output[1]} -o {output[0]} {params.cyto_arg} &> {log}"
