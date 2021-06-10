import glob
import os
import sys

sys.path.insert(0, workflow.basedir)
from constants import VERSION

print(f"Ploidetect-pipeline {VERSION}")


## Load default config values
configfile: os.path.join(workflow.basedir, "resources/config/default_run_params.yaml")
configfile: os.path.join(workflow.basedir, "resources/config/default_case.yaml")
configfile: os.path.join(workflow.basedir, "resources/config/genome_ref.yaml")


MEM_PER_CPU = 7900

chromosomes = config["chromosomes"]
output_dir = config["output_dir"]

scripts_dir = os.path.join(workflow.basedir, "scripts")
array_positions = (
    config["array_positions"][config["genome_name"]]
    if os.path.exists(config["array_positions"][config["genome_name"]])
    else os.path.join(
        workflow.basedir, config["array_positions"][config["genome_name"]]
    )
)


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


def devtools_install():
    if config["ploidetect_local_clone"] and config["ploidetect_local_clone"] != "None":
        install_path = config["ploidetect_local_clone"].format(**config)
        devtools_cmd = f"devtools::install_local('{install_path}', force = TRUE)"
    else:
        ver = config["ploidetect_ver"]
        devtools_cmd = f"devtools::install_github('lculibrk/Ploidetect', ref = '{ver}')"
    return f'"{devtools_cmd}"'


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
    shell:
        "export LC_ALL=en_US.UTF-8; "
        " Rscript -e {params} "
        " && echo {params} > {output} && date >> {output}"


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
    log:
        "{output_dir}/logs/germline_cov.{case}.{normal}.{chr}.log",
    shell:
        "samtools depth -r{wildcards.chr} -Q 21 {input.bam} 2>> {log}"
        " | awk -v FS='\t' -v OFS='\t' 'NR > 1{{print $1, $2, $2+1, $3}}'"
        " | python3 {params.scripts_dir}/make_windows.py - {params.threshold} 2>> {log}"
        " | bedtools sort -i stdin > {output}  2>> {log}"
        " && ls -l {output} >> {log}"


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
    """Compute depth of tumor and normal reads in previously created bins"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        nombam=lambda w: config["bams"][w.case]["normal"][w.normal],
        window=rules.splitwindowfile.output,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/tumour/{chr}.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/genomecovsomatic.{case}.{somatic}_{normal}.{chr}.log",
    shell:
        "bedtools multicov -bams {input.sombam} {input.nombam} -q 20 -bed {input.window}  > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule mergesomatic:
    """Merge output of genomecovsomatic to a singular file"""
    input:
        expand(
            "{{output_dir}}/scratch/{{case}}/{{somatic}}_{{normal}}/tumour/{chr}.bed",
            chr=chromosomes,
        ),
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/tumour.bed"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/mergesomatic.{case}.{somatic}_{normal}.log",
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule compute_loh:
    """Variant call the germline, filter for heterozygous snps and count alleles in somatic"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        normbam=lambda w: config["bams"][w.case]["normal"][w.normal],
    output:
        folder=directory("{output_dir}/scratch/{case}/{somatic}_{normal}/loh_tmp"),
        loh=temp("{output_dir}/scratch/{case}/{somatic}_{normal}/loh_tmp/loh_raw.txt"),
    params:
        genome=config["genome"][config["genome_name"]],
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
        " {params.genome} {params.array_positions}"
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
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/gc.bed"),
    params:
        genome=config["genome"][config["genome_name"]],
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
        "bedtools nuc -fi {params.genome} -bed {input} | cut -f1,2,3,5 | tail -n +2 > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule mergedbed:
    """Merge all the data into a single file"""
    input:
        gc=rules.getgc.output,
        tumour=rules.mergesomatic.output,
        normal=rules.merge_germline.output,
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
        "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc}) > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        "{output_dir}/scratch/{case}/{somatic}_{normal}/merged.bed",
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
        "Rscript {params.scripts_dir}/prep_ploidetect2.R -i {input} -o {output} &> {log}"


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
    log:
        "{output_dir}/logs/ploidetect.{case}.{somatic}_{normal}.log",
    shell:
        "Rscript {params.scripts_dir}/run_ploidetect2.R "
        " -i {input.preseg} "
        " -o {output.models} -p {output.plots} -r {output.meta}"
        " &> {log}"


rule ploidetect_copynumber:
    """Performs CNV calling using the tumor purity and ploidy estimated by Ploidetect"""
    input:
        models="{output_dir}/{case}/{somatic}_{normal}/models.txt",
        segs="{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
        ploidetect_plots="{output_dir}/{case}/{somatic}_{normal}/plots.pdf",
    output:
        cna="{output_dir}/{case}/{somatic}_{normal}/cna.txt",
        cna_plots="{output_dir}/{case}/{somatic}_{normal}/cna_plots.pdf",
        cna_cond="{output_dir}/{case}/{somatic}_{normal}/cna_condensed.txt",
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    params:
        scripts_dir=scripts_dir,
    log:
        "{output_dir}/logs/ploidetect_copynumber.{case}.{somatic}_{normal}.log",
    shell:
        "Rscript {params.scripts_dir}/ploidetect_copynumber.R"
        " -i {input.segs} -m {input.models}"
        " -p {output.cna_plots} -o {output.cna}"
        " &> {log}"
