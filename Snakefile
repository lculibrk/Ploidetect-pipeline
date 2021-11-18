import glob
import os
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()
import requests
HTTP = HTTPRemoteProvider()

sys.path.insert(0, workflow.basedir)
from constants import VERSION
from constants import WorkflowSetupError

print(f"Ploidetect-pipeline {VERSION}")


## Load config values
configfile: os.path.join(workflow.basedir, "config/defaults.yaml")
configfile: os.path.join(workflow.basedir, "config/samples.yaml")
configfile: os.path.join(workflow.basedir, "config/parameters.yaml")

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

scripts_dir = os.path.join(workflow.basedir, "scripts")
array_positions = (
    config["array_positions"][config["genome_name"]]
    if os.path.exists(config["array_positions"][config["genome_name"]])
    else os.path.join(
        workflow.basedir, config["array_positions"][config["genome_name"]]
    )
)

with open(array_positions) as f:
    snp_positions = f.readlines()
    snp_positions = [snp.strip().split("\t") for snp in snp_positions]
    snp_chrs = [line[0] for line in snp_positions]
    snp_chrs = list(set(snp_chrs))
    for chromosome in chromosomes:
        if str(chromosome) not in snp_chrs:
            raise(WorkflowSetupError(f"Chromosome {chromosome} specified in config was not found in provided snp position data"))
    lens = [len(line) for line in snp_positions]
    print(lens[0])
    if any(leng > 2 for leng in lens):
        raise(WorkflowSetupError("More than two columns detected in snp position data - file must contain only chromosome and position columns"))
    



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
    somatic_paths = bams_dict[sample]["somatic"].values()
    normals = bams_dict[sample]["normal"].keys()
    normal_paths = bams_dict[sample]["normal"].values()
    ## Check that the bams are findable
    concat_bams = list(somatic_paths) + list(normal_paths)
#    for bam in concat_bams:
#        if not os.path.exists(bam):
#            #raise WorkflowSetupError(f"Input file {bam} could not be found. Ensure that the file is spelled correctly, and that you've correctly bound the directory if using singularity")
    combinations = expand("{somatic}_{normal}", somatic=somatics, normal=normals)
    #outs = [os.path.join(output_dir, sample, comb, "cna.txt") for comb in combinations]
    outs = [os.path.join(output_dir, "scratch", sample, comb, "merged.bed") for comb in combinations]
    output_list.extend(outs)

print(f"Final outputs: {output_list}")

rule all:
    input:
        output_list


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
        install_cmd = devtools_install(),
    log:
        expand("{output_dir}/logs/ploidetect_install.log", output_dir = output_dir)
    benchmark:
        "benchmark/pdt_install.txt"
    shell:
        "export LC_ALL=en_US.UTF-8; "
        " Rscript -e {params.install_cmd} > {log}"
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
    benchmark:
        "benchmark/cytobands.txt"
    shell:
        "gunzip -c {input} | sed 's/chr//g' > {output}"

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
        "docker://lculibrk/ploidetect"
    shell:
        "gunzip -c {input} > {output[0]}; samtools faidx {output[0]}"

rule copy_ncram:
    """Copy over germline cram(i)"""
    input:
        cram=lambda w: GS.remote(config["bams"][w.case]["normal"][w.lib]),
        crai=lambda w: GS.remote(config["bams"][w.case]["normal"][w.lib] + ".crai")
    output:
        cram=temp("{output_dir}/scratch/{case}/{lib}/normal.cram"),
        crai=temp("{output_dir}/scratch/{case}/{lib}/normal.cram.crai")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    benchmark:
        "{output_dir}/benchmark/{case}/{lib}/copy_normal.txt"
    shell:
        "cp {input.cram} {output.cram}; cp {input.crai} {output.crai}"


rule copy_cram:
    """Copy over somatic cram(i)"""
    input:
        cram=lambda w: GS.remote(config["bams"][w.case]["somatic"][w.lib]),
        crai=lambda w: GS.remote(config["bams"][w.case]["somatic"][w.lib] + ".crai")
    output:
        cram=temp("{output_dir}/scratch/{case}/{lib}/somatic.cram"),
        crai=temp("{output_dir}/scratch/{case}/{lib}/somatic.cram.crai")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    benchmark:
        "{output_dir}/benchmark/{case}/{lib}/copy_tumour.txt"
    shell:
        "cp {input.cram} {output.cram}; cp {input.crai} {output.crai}"


rule germline_cov:
    """Compute per-base depth in germline bam"""
    input:
        bam="{output_dir}/scratch/{case}/{normal}/normal.cram",
        bami="{output_dir}/scratch/{case}/{normal}/normal.cram.crai",
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
        qual = config[config["sequence_type"]]["qual"],
        maxd = config[config["sequence_type"]]["maxd"]
    log:
        "{output_dir}/logs/germline_cov.{case}.{normal}.{chr}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{normal}/germline_cov{chr}.txt"
    shell:
        "samtools depth -r{wildcards.chr} -Q{params.qual} -m {params.maxd} {input.bam} > {output} 2>> {log}"

rule tumour_cov:
    """Compute per-base depth in tumour bam"""
    input:
        bam="{output_dir}/scratch/{case}/{tumour}/somatic.cram",
        bami="{output_dir}/scratch/{case}/{tumour}/somatic.cram.crai",
    output:
        temp("{output_dir}/scratch/{case}/{tumour}/tumour/{chr}.bed"),
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
        qual = config[config["sequence_type"]]["qual"],
        maxd = config[config["sequence_type"]]["maxd"]
    log:
        "{output_dir}/logs/germline_cov.{case}.{tumour}.{chr}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{tumour}/germline_cov{chr}.txt"
    shell:
        "samtools depth -r{wildcards.chr} -Q{params.qual} -m {params.maxd} {input.bam} > {output} 2>> {log}"

rule concat_tumour:
    input:
        expand("{{output_dir}}/scratch/{{case}}/{{lib}}/tumour/{chr}.bed", chr = chromosomes)
    output:
        temp("{output_dir}/scratch/{case}/{lib}/tumour.bed")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/concat_depths.{case}.{lib}_tumour.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{lib}/concat_tumour.txt"
    shell:
        "cat {input} > {output}"

rule concat_normal:
    input:
        expand("{{output_dir}}/scratch/{{case}}/{{lib}}/normal/{chr}.bed", chr = chromosomes)
    output:
        temp("{output_dir}/scratch/{case}/{lib}/normal.bed")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/concat_normal.{case}.{lib}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{lib}/concat_normal.txt"
    shell:
        "cat {input} > {output}"
        
rule make_bins:
    input:
        "{output_dir}/scratch/{case}/{normal}/normal/normal.bed"
    output:
        temp("{output_dir}/scratch/{case}/{normal}/normal/windows_unsort.txt")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/make_bins.{case}.{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{normal}/make_bins.txt"
    shell:
        "python3 {params.scripts_dir}/make_windows.py {input} {params.threshold} > {output} 2>> {log}"

rule sort_bins:
    input:
        "{output_dir}/scratch/{case}/{normal}/normal/windows_unsort.txt"
    output:
        temp("{output_dir}/scratch/{case}/{normal}/normal/windows.txt")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/sort_bins.{case}.{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{normal}/sort_bins.txt"
    shell:
        "sort -k1,1 -k2,2n -S 300M -T tmp/ {input} > {output} 2>> {log}"
        
rule genomecovsomatic:
    input:
        depth="{output_dir}/scratch/{case}/{somatic}/tumour.bed",
        window=rules.sort_bins.output
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/{somatic}.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/genomecovsomatic.txt"
    shell:
        "python3 scripts/summarize_counts.py {input.depth} {input.window} > {output}"
            
rule genomecovnormal:
    input:
        depth="{output_dir}/scratch/{case}/{normal}/normal.bed",
        window=rules.sort_bins.output
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/{normal}.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/genomecovnormal.txt"
    shell:
        "python3 scripts/summarize_counts.py {input.depth} {input.window} > {output}"
    

rule split_positions:
    """Split the array positions file by chromosome for parallel processing"""
    input:
        array_positions
    output:
        temp("{output_dir}/scratch/split_array/{chr}.txt")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/split_array_{chr}.log",
    benchmark:
        "{output_dir}/benchmark/splitpositions{chr}.txt"
    shell:
        "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"
    
rule generate_bafs:
    input:
        normbam="{output_dir}/scratch/{case}/{normal}/normal.cram",
        nori="{output_dir}/scratch/{case}/{normal}/normal.cram.crai",
        sombam="{output_dir}/scratch/{case}/{somatic}/somatic.cram",
        somi="{output_dir}/scratch/{case}/{somatic}/somatic.cram.crai",
        genome=genome_path,
        array="{output_dir}/scratch/split_array/{chr}.txt"
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/bafs/{chr}.txt")
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
        array_positions=array_positions,
    log:
        "{output_dir}/logs/generate_bafs.{case}.{somatic}_{normal}_{chr}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/generate_bafs_{chr}.txt"
    shell:
        "python {params.scripts_dir}/light_varcaller.py {input.normbam} {input.sombam} {input.array} {input.genome} > {output}"
    
rule concat_bafs:
    input:
        expand("{{output_dir}}/scratch/{{case}}/{{somatic}}_{{normal}}/bafs/{chr}.txt", chr = chromosomes)
    output:
        "{output_dir}/scratch/{case}/{somatic}_{normal}/bafs.txt"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    container:
        "docker://lculibrk/ploidetect"
    params:
        scripts_dir=scripts_dir,
        array_positions=array_positions,
    log:
        "{output_dir}/logs/concat_bafs.{case}.{somatic}_{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/concat_bafs.txt"
    shell:
        "sort -k1,1 -k2,2n -S 300M -T tmp/ {input} > {output}"

rule getgc:
    """Get GC content for each bin"""
    input:
        window=rules.sort_bins.output,
        genome = genome_path,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}_{normal}/gc.bed"),
    params:
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/getgc.{case}.{somatic}_{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/getgc.txt"
    shell:
        "python3 scripts/stream_gc.py {input.genome} {input.window} > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule mergedbed:
    """Merge all the data into a single file"""
    input:
        gc=rules.getgc.output,
        tumour="{output_dir}/scratch/{case}/{somatic}_{normal}/tumour.bed",
        normal="{output_dir}/scratch/{case}/{somatic}_{normal}/normal.bed",
        loh=rules.concat_bafs.output,
    output:
        "{output_dir}/scratch/{case}/{somatic}_{normal}/merged.bed",
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/mergedbed.{case}.{somatic}_{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/mergedbed.txt"
    shell:
        "paste {input.tumour} <(cut -f4 {input.normal}) <(cut -f4 {input.loh}) <(cut -f4 {input.gc})"
        "| sed 's/chr//g' > {output}" ## Cuts out any "chr" if using hg38
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule preseg:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        "{output_dir}/scratch/{case}/{somatic}_{normal}/merged.bed",
        rules.ploidetect_install.output if workflow.use_conda and not workflow.use_singularity else __file__,
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
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/preseg.txt"
    shell:
        "Rscript {scripts_dir}/prep_ploidetect2.R -i {input[0]} -c {input.cytos} -o {output}"

        
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
        )
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
        cyto_arg = cyto_arg,
    log:
        "{output_dir}/logs/ploidetect.{case}.{somatic}_{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/ploidetect.txt"
    shell:
        "Rscript {params.scripts_dir}/run_ploidetect2.R "
        " -i {input.preseg} "
        " -o {output.models} -p {output.plots} -r {output.meta} {params.cyto_arg}"
        " &> {log}"


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
        "docker://lculibrk/ploidetect"
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    params:
        scripts_dir=scripts_dir,
        cyto_arg = cyto_arg,
    log:
        "{output_dir}/logs/ploidetect_copynumber.{case}.{somatic}_{normal}.log",
    benchmark:
        "{output_dir}/benchmark/{case}/{somatic}_{normal}/cnv.txt"
    shell:
        "Rscript {params.scripts_dir}/ploidetect_copynumber.R"
        " -i {input.rds} -m {input.models}"
        " -p {output.cna_plots} -o {output.cna} {params.cyto_arg}"
        " &> {log}"
