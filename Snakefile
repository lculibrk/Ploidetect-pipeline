import glob
import os
import sys

sys.path.insert(0, workflow.basedir)
from constants import VERSION

print(f"Ploidetect-pipeline {VERSION}")


## Load config values
configfile: os.path.join(workflow.basedir, "CONFIG.txt")



MEM_PER_CPU = 7900

chromosomes = config["chromosomes"]
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

def nanopore_handling():
    if config["sequence_type"] == "ont":
        return {"maxd":500, "qual":20}
    else:
        return {"maxd":0,   "qual":50}


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
    params: 
        qual = nanopore_handling()["qual"],
        maxd = nanopore_handling()["maxd"]
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "samtools depth -Q{params.qual} -m {params.maxd} -r{wildcards.chr} {input[0]}"
        " | awk -v FS='\\t' -v OFS='\\t' 'NR > 1{{print $1, $2, $2+1, $3}}'"
        " | python3 scripts/make_windows.py - 100000  > {output}"



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
        "docker://lculibrk/ploidetect"
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
        "docker://lculibrk/ploidetect"
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
        "docker://lculibrk/ploidetect"
    shell:
        "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"

rule genomecovsomatic:
    input:
        lambda w: config["bams"][w.case]["somatic"][w.somatic],
        window=rules.splitwindowfile.output
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/tumour/{chr}.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU
    params: 
        qual = nanopore_handling()["qual"],
        maxd = nanopore_handling()["maxd"]
    shell:
        "samtools depth -Q {params.qual} -m {params.maxd} -r {wildcards.chr} -a {input[0]} "
        " | awk -v FS='\\t' -v OFS='\\t' \'{{print $1, $2, $2 + 1, $3}}\'"
        " | sort -k1,1 -k2,2n "
        " | bedtools map -b stdin -a {input.window} -c 4 -o mean > {output}"

rule genomecovgermline:
    input:
        lambda w: config["bams"][w.case]["normal"][w.normal],
        window=rules.splitwindowfile.output
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/normal/{chr}.bed")
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    params: 
        qual = nanopore_handling()["qual"],
        maxd = nanopore_handling()["maxd"]
    shell:
        "samtools depth -Q {params.qual} -m {params.maxd} -r {wildcards.chr} -a {input[0]}"
        " | awk -v FS='\\t' -v OFS='\\t' \'{{print $1, $2, $2 + 1, $3}}\'"
        " | sort -k1,1 -k2,2n"
        " | bedtools map -b stdin -a {input.window} -c 4 -o mean > {output}"


rule merge_split_tumour:
    input:
        expand("{{temp_dir}}/{{case}}/{{somatic}}_{{normal}}/tumour/{chr}.bed", chr = chromosomes)
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/tumour.bed")
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    shell:
        "cat {input}"
        " | sort -k1,1 -k2,2n > {output}"


rule merge_split_normal:
    input:
        expand("{{temp_dir}}/{{case}}/{{somatic}}_{{normal}}/normal/{chr}.bed", chr = chromosomes)
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/normal.bed")
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    shell:
        "cat {input}"
        " | sort -k1,1 -k2,2n > {output}"


rule compute_loh:
    """Variant call the germline, filter for heterozygous snps and count alleles in somatic"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        normbam=lambda w: config["bams"][w.case]["normal"][w.normal],
    output:
        temp(directory("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp/")),
        temp("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp/loh_raw.txt"),
    params:
        genome=config["genome"][config["genome_name"]],
        array_positions={array_positions},
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "bash {scripts_dir}/get_allele_freqs.bash {input.normbam} {input.sombam}"
        " {params.genome} {params.array_positions}"
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
        "docker://lculibrk/ploidetect"
    shell:
        "awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}' {input.loh}"
        " | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}'"
        " | bedtools sort -i stdin"
        " | Rscript {scripts_dir}/merge_loh.R -l STDIN -w {input.window} -o {output}"


rule getgc:
    """Get GC content for each bin"""
    input:
        window=rules.makewindowfile.output,
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/gc.bed"),
    params:
        genome=config["genome"][config["genome_name"]],
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "bedtools nuc -fi {params.genome} -bed {input} | cut -f1,2,3,5 | tail -n +2 > {output}"


rule mergedbed:
    """Merge all the data into a single file"""
    input:
        gc=rules.getgc.output,
        tumour="{temp_dir}/{case}/{somatic}_{normal}/tumour.bed",
        normal="{temp_dir}/{case}/{somatic}_{normal}/normal.bed",
        loh=rules.process_loh.output,
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/merged.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "paste {input.tumour} <(cut -f4 {input.normal}) <(cut -f4 {input.loh}) <(cut -f4 {input.gc}) > {output}"


rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        temp_dir + "{case}/{somatic}_{normal}/merged.bed",
    output:
        "{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "Rscript {scripts_dir}/prep_ploidetect2.R -i {input} -o {output}"


rule ploidetect:
    """Runs Ploidetect"""
    input:
        rules.preseg.output,
        rules.ploidetect_install.output if not workflow.use_singularity and "install_ploidetect" in config.keys() and config[
            "install_ploidetect"
        ] else __file__,
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
    shell:
        "Rscript {scripts_dir}/run_ploidetect2.R "
        " -i {input[0]} "
        " -o {output.models} -p {output.plots} -r {output.meta}"


rule ploidetect_copynumber:
    """Performs CNV calling using the tumor purity and ploidy estimated by Ploidetect"""
    input:
        "{output_dir}/{case}/{somatic}_{normal}/models.txt",
        "{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
        "{output_dir}/{case}/{somatic}_{normal}/plots.pdf",
    output:
        "{output_dir}/{case}/{somatic}_{normal}/cna.txt",
        "{output_dir}/{case}/{somatic}_{normal}/cna_plots.pdf",
        "{output_dir}/{case}/{somatic}_{normal}/cna_condensed.txt",
    conda:
        "conda_configs/r.yaml"
    log:
        "{output_dir}/{case}/{somatic}_{normal}/cna_log.txt"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    shell:
        "Rscript {scripts_dir}/ploidetect_copynumber.R -i {input[1]} -m {input[0]} -p {output[1]} -o {output[0]} &> {log}"
