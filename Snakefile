import glob
import os

## Versioning
__version__ = "v0.0.3"  # merge - GSC and multi-sample
print(f"Ploidetect-pipeline {__version__}")

# Use --configfile to specify - test default shown below
if "chromosomes" not in config.keys():
    configfile: os.path.join(workflow.basedir, "CONFIG.txt")

## Load config values
chromosomes=config["chromosomes"]
output_dir = config["output_dir"]
temp_dir = config["temp_dir"] if config["temp_dir"] else f"{output_dir}/temp/"
if temp_dir[-1] != "/":
    temp_dir += "/"  # Prevents strange case wild-card error.

## Parse sample information
bams_dict = config['bams']
sample_ids = bams_dict.keys()

## Here comes the for loop
## Initialize output list, which will be fed to rule all
output_list = []
## Loop over the bams and construct lib comparisons
for sample in sample_ids:
    somatics = bams_dict[sample]["somatic"].keys()
    normals = bams_dict[sample]["normal"].keys()
    combinations = expand("{somatic}_{normal}", somatic = somatics, normal = normals)
    outs = [os.path.join(output_dir, sample, comb, "cna.txt") for comb in combinations]
    output_list.extend(outs)



scripts_dir = os.path.join(workflow.basedir, "scripts")
array_positions = config["array_positions"][config["genome_name"]] if os.path.exists(config["array_positions"][config["genome_name"]]) else os.path.join(workflow.basedir, config["array_positions"][config["genome_name"]])

print(output_list)

rule all:
    input:
        output_list

def check_docker():
    """Ploidetect installed filepath.
    Filename to create on a successful install or check for successful installation.
    Checks the config file for docker options.
    """
    if config["use-docker"] == 1:
    # /dev/null should be present in basically every 'nix system
    # This ensures that install_ploidetect isn't run
        file_to_make = "/dev/null"
    else:
    # if not using docker, should install Ploidetect
        file_to_make = os.path.join(workflow.basedir, "conda_configs/ploidetect_installed.txt")
    return(file_to_make)

def devtools_install():
    if config["ploidetect_local_clone"] and config["ploidetect_local_clone"] != "None":
        install_path = config["ploidetect_local_clone"].format(**config)
        devtools_cmd = "\"devtools::install_local('" + install_path + "')\""
    else:
        devtools_cmd = "\"devtools::install_github('lculibrk/Ploidetect', "
        devtools_cmd += "ref = '" + config["ploidetect_github_version"] + "')\""
    return(devtools_cmd)

rule install_ploidetect:
    """Install Ploidetect R script into environment"""
    output:
        expand("{install_dir}/conda_configs/ploidetect_installed.txt", install_dir=workflow.basedir)
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/r.yaml"
    params:
        devtools_install()
    shell:
        "export LC_ALL=en_US.UTF-8; "
        " Rscript -e {params} "
        " && echo {params} > {output} && date >> {output}"


rule germline_cov:
    """Compute per-base depth in germline bam, convert to .bed format and pile up into equal-coverage bins"""
    input:
        bam = lambda w: config["bams"][w.case]["normal"][w.lib],
    output:
        temp("{temp_dir}/{case}/{lib}/normal/{chr}.bed")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "samtools depth -r{wildcards.chr} -Q 21 {input.bam}"
        " | awk -v FS='\t' -v OFS='\t' 'NR > 1{{print $1, $2, $2+1, $3}}'"
        " | python3 {scripts_dir}/make_windows.py - 100000"
        " | bedtools sort -i stdin > {output}"

rule merge_germline:
    """Merge multi-chromosome output from germline_cov into single file"""
    input:
        expand("{{temp_dir}}/{{case}}/{{normal}}/normal/{chr}.bed", chr=chromosomes)
    output:
        temp("{temp_dir}/{case}/{normal}/germline.bed")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"

rule makewindowfile:
    """Remove germline depth column from file to obtain bins"""
    input:
        rules.merge_germline.output
    output:
        temp("{temp_dir}/{case}/{normal}/windows.txt")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "cut -f 1,2,3 < {input} | bedtools sort -i stdin > {output}"

rule splitwindowfile:
    """Split bins into each chromosome for parallel computing of depth in somatic"""
    input:
        rules.makewindowfile.output
    output:
        temp("{temp_dir}/{case}/{normal}/windows/{chr}.txt")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"

rule genomecovsomatic:
    """Compute depth of tumor and normal reads in previously created bins"""
    input:
        sombam= lambda w: config["bams"][w.case]["somatic"][w.somatic],
        nombam= lambda w: config["bams"][w.case]["normal"][w.normal],
        window=rules.splitwindowfile.output
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/tumour/{chr}.bed")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "bedtools multicov -bams {input.sombam} {input.nombam} -q 20 -bed {input.window}  > {output}"

rule mergesomatic:
    """Merge output of genomecovsomatic to a singular file"""
    input:
        expand("{{temp_dir}}/{{case}}/{{somatic}}_{{normal}}/tumour/{chr}.bed", chr=chromosomes)
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/tumour.bed")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"


rule compute_loh:
    """Variant call the germline, filter for heterozygous snps and count alleles in somatic"""
    input:
        sombam= lambda w: config["bams"][w.case]["somatic"][w.somatic],
        normbam= lambda w: config["bams"][w.case]["normal"][w.normal],
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/loh_tmp/loh_raw.txt")
    params:
        genome = config["genome"][config["genome_name"]],
	    array_positions = {array_positions}
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "bash {scripts_dir}/get_allele_freqs.bash {input.normbam} {input.sombam}"
        " {params.genome} {params.array_positions} {temp_dir}/loh_tmp/"

rule process_loh:
    """Convert allele counts to beta-allele frequencies and merge for each bin"""
    input:
        loh=rules.compute_loh.output,
        window=rules.makewindowfile.output
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/loh.bed")
    resources: cpus=1, mem_mb=7900
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
        window=rules.makewindowfile.output
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/gc.bed")
    params:
        genome=config["genome"][config["genome_name"]]
    resources: cpus=1, mem_mb=7900
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
	    loh=rules.process_loh.output
    output:
        temp("{temp_dir}/{case}/{somatic}_{normal}/merged.bed")
    resources: cpus=1, mem_mb=7900
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc}) > {output}"

rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        temp_dir + "{case}/{somatic}_{normal}/merged.bed"
    output:
        "{output_dir}/{case}/{somatic}_{normal}/segmented.RDS"
    resources: cpus=24, mem_mb=189600
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
	    check_docker()
    output:
        plots="{output_dir}/{case}/{somatic}_{normal}/plots.pdf",
        models="{output_dir}/{case}/{somatic}_{normal}/models.txt",
        meta="{output_dir}/{case}/{somatic}_{normal}/meta.RDS"
    conda:
        "conda_configs/r.yaml"
    resources: cpus=24, mem_mb=189600
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
        "{output_dir}/{case}/{somatic}_{normal}/plots.pdf"
    output:
        "{output_dir}/{case}/{somatic}_{normal}/cna.txt",
        "{output_dir}/{case}/{somatic}_{normal}/cna_plots.pdf"
    conda:
        "conda_configs/r.yaml"
    log:
        "{output_dir}/{case}/{somatic}_{normal}/cna_log.txt"
    resources: cpus=24, mem_mb=189600
    shell:
        "Rscript {scripts_dir}/ploidetect_copynumber.R -i {input[1]} -m {input[0]} -p {output[1]} -o {output[0]} &> {log}"
