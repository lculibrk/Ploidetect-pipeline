rule mkgenome_ss:
    """Make bedtools makewindows input"""
    input:
        genome = genome_path + ".fai"
    output:
        temp("{output_dir}/scratch/genome.bins")
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "cut -f 1,2 {input} > {output}"

rule filter_windows:
    """Filters out non-standard chromosomes from windows file"""
    input:
        rules.mkgenome_ss.output,
    output:
        temp("{output_dir}/scratch/genome.filt")
    run:
        ps = "{{print $0}}"
        for chrom in chromosomes:
            print("awk '$1 == \"" + chrom + "\" " + ps + "' " + input[0] + " >> " + output[0])
            shell("awk '$1 == \"" + chrom + "\" " + ps + "' " + input[0] + " >> " + output[0])

rule windows_ss:
    """Create bins for single-sample data"""
    input:
        genome = "{output_dir}/scratch/genome.filt"
    output:
        temp("{output_dir}/scratch/{case}/{somatic}/windows.bed"),
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "bedtools makewindows -g {input.genome} -w 5000 > {output}"



rule splitwindowfile_ss:
    """Split bins by chromosome for paralellization"""
    input:
        rules.windows_ss.output,
    output:
        temp("{output_dir}/scratch/{case}/{somatic}/windows/{chr}.txt"),
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/splitwindowfile.{case}.{somatic}.{chr}.log",
    shell:
        "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"

rule genomecovsomatic_ss:
    input:
        lambda w: config["bams"][w.case]["somatic"][w.somatic],
        window=rules.splitwindowfile_ss.output
    output:
        temp("{output_dir}/scratch/{case}/{somatic}/tumour/{chr}.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU
    params: 
        qual = config[config["sequence_type"]]["qual"],
        maxd = config[config["sequence_type"]]["maxd"]
    wildcard_constraints:
        soamtic="[^_]*"
    shell:
        "samtools depth -Q {params.qual} -m {params.maxd} -r {wildcards.chr} -a {input[0]} "
        " | awk -v FS='\\t' -v OFS='\\t' \'{{print $1, $2, $2 + 1, $3}}\'"
        " | sort -k1,1 -k2,2n "
        " | bedtools map -b stdin -a {input.window} -c 4 -o mean > {output}"

rule merge_split_tumour_ss:
    input:
        expand(
            "{{output_dir}}/scratch/{{case}}/{{somatic}}/tumour/{chr}.bed",
            chr=chromosomes,
        )
    output:
        temp("{output_dir}/scratch/{case}/{somatic}/tumour.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    log:
        "{output_dir}/logs/mergesomatic.{case}.{somatic}.log"
    shell:
        "cat {input} | bedtools sort -i stdin > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"

rule loh_ss:
    """Compute allele frequencies for single-sample data"""
    input:
        sombam=lambda w: config["bams"][w.case]["somatic"][w.somatic],
        window=rules.windows_ss.output,
        genome = genome_path,
    output:
        "{output_dir}/scratch/{case}/{somatic}/loh.bed",
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    params:
        array_positions={array_positions},
        scripts_dir=scripts_dir,
    log:
        "{output_dir}/logs/process_loh.{case}.{somatic}.log",
    shell:
        "samtools mpileup {input.sombam} -l {params.array_positions} -f {input.genome} -B"
        " | awk -f {params.scripts_dir}/parse_pileup.awk"
        " | awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}'"
        " | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}'"
        " | bedtools sort -i stdin"
        " | Rscript {params.scripts_dir}/merge_loh.R -l STDIN -w {input.window} -o {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"


rule getgc_ss:
    """Get GC content for each bin"""
    input:
        window=rules.windows_ss.output,
        genome = genome_path,
    output:
        "{output_dir}/scratch/{case}/{somatic}/gc.bed",
    params:
    resources:
        cpus=1,
        mem_mb=MEM_PER_CPU,
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/getgc.{case}.{somatic}.log",
    shell:
        "bedtools nuc -fi {input.genome} -bed {input.window} | cut -f1,2,3,5 | tail -n +2 > {output}"
        " 2> {log}"
        " && ls -l {output} >> {log}"

rule merge_single:
    """Merge single-sample data into a single file"""
    input:
        gc="{output_dir}/scratch/{case}/{somatic}/gc.bed",
        tumour="{output_dir}/scratch/{case}/{somatic}/tumour.bed",
        loh="{output_dir}/scratch/{case}/{somatic}/loh.bed",
    output:
        temp("{output_dir}/scratch/{case}/{somatic}/merged.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/mergedbed.{case}.{somatic}.log",
    shell:
        "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc})"
        "| sed 's/chr//g' > {output}" ## Cuts out any "chr" if using hg38
        " 2> {log}"
        " && ls -l {output} >> {log}"

rule preseg_ss:
    """Presegment and prepare data for input into Ploidetect"""
    input:
        "{output_dir}/scratch/{case}/{somatic}/merged.bed",
        rules.ploidetect_install.output if workflow.use_conda and not workflow.use_singularity else __file__,
        cytos=cyto_path,
    output:
        "{output_dir}/{case}/{somatic}/segmented.RDS",
    resources:
        cpus=24,
        mem_mb=24 * MEM_PER_CPU,
    conda:
        "conda_configs/r.yaml"
    container:
        "docker://lculibrk/ploidetect:v1.4.0"
    params:
        scripts_dir=scripts_dir,
    log:
        "{output_dir}/logs/preseg.{case}.{somatic}.log",
    shell:
        "Rscript {scripts_dir}/prep_ploidetect2.R -i {input[0]} -c {input.cytos} -o {output} -s TRUE"
