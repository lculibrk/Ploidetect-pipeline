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
                "{output_dir}/{sample}/{somatic}_{normal}/cna.txt",
                output_dir=output_dir,
                sample=sample,
                somatic=config["bams"][sample]["somatic"].keys(),
                normal=config["bams"][sample]["normal"].keys(),
            )
            for sample in config["bams"].keys()
        ],


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
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
    shell:
        "samtools depth -r{wildcards.chr} -Q 21 {input.bam}"
        " | awk -v FS='\t' -v OFS='\t' 'NR > 1{{print $1, $2, $2+1, $3}}'"
        " | python3 {params.scripts_dir}/make_windows.py - 100000"
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
        "docker://lculibrk/ploidetect"
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
        "docker://lculibrk/ploidetect"
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"


rule compute_loh:
    """Variant call the germline, filter for heterozygous snps and count alleles in somatic"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        normbam=lambda w: config["bams"][w.case]["normal"][w.normal],
    output:
        folder=directory("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp"),
        loh=temp("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp/loh_raw.txt"),
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
    shell:
        "bash {params.scripts_dir}/get_allele_freqs.bash {input.normbam} {input.sombam}"
        " {params.genome} {params.array_positions}"
        " {output.folder}"


rule process_loh:
    """Convert allele counts to beta-allele frequencies and merge for each bin"""
    input:
        loh=rules.compute_loh.output.loh,
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
    params:
        scripts_dir=scripts_dir,
    shell:
        "awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}' {input.loh}"
        " | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}'"
        " | bedtools sort -i stdin"
        " | Rscript {params.scripts_dir}/merge_loh.R -l STDIN -w {input.window} -o {output}"


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
        "docker://lculibrk/ploidetect"
    shell:
        "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc}) > {output}"


rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    # Use explict params to connect temp_dir and output_dir wildcards
    input:
        #       "{temp_dir}/{case}/{somatic}_{normal}/merged.bed",
        [
            expand(
                "{temp_dir}/{sample}/{somatic}_{normal}/merged.bed",
                output_dir=output_dir,
                sample=sample,
                somatic=config["bams"][sample]["somatic"].keys(),
                normal=config["bams"][sample]["normal"].keys(),
                temp_dir=temp_dir,
            )
            for sample in config["bams"].keys()
        ],
    output:
        #       "{output_dir}/{case}/{somatic}_{normal}/segmented.RDS",
        [
            expand(
                "{output_dir}/{sample}/{somatic}_{normal}/segmented.RDS",
                output_dir=output_dir,
                sample=sample,
                somatic=config["bams"][sample]["somatic"].keys(),
                normal=config["bams"][sample]["normal"].keys(),
                temp_dir=temp_dir,
            )
            for sample in config["bams"].keys()
        ],
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect"
    script:
        "scripts/prep_ploidetect2.R -i {input} -o {output}"


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
    script:
        "scripts/run_ploidetect2.R "
        " -i {input.preseg} "
        " -o {output.models} -p {output.plots} -r {output.meta}"


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
    log:
        "{output_dir}/{case}/{somatic}_{normal}/cna_log.txt",
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    script:
        "scripts/ploidetect_copynumber.R -i {input.segs} -m {input.models} -p {output.cna_plots} -o {output.cna} &> {log}"
