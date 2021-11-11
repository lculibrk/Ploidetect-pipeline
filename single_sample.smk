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

rule windows_ss:
    """Create bins for single-sample data"""
    input:
        genome = "{output_dir}/scratch/genome.bins"
    output:
        temp("{output_dir}/scratch/{case}/{somatic}/windows.bed"),
    container:
        "docker://lculibrk/ploidetect"
    shell:
        "bedtools makewindows -g {input.genome} -w 100000 > {output}"

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
        temp("{output_dir}/scratch/{case}/{somatic}/loh.bed"),
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
        " | bedtools sort -i stdin"
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
        temp("{output_dir}/scratch/{case}/{somatic}/gc.bed"),
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
        gc="{output_dir}/scratch/{case}/{tumour}/gc.bed",
        tumour="{output_dir}/scratch/{case}/{tumour}/tumour.bed",
        loh="{output_dir}/scratch/{case}/{tumour}/loh.bed",
    output:
        temp("{output_dir}/scratch/{case}/{tumour}/merged.bed"),
    conda:
        "conda_configs/sequence_processing.yaml"
    container:
        "docker://lculibrk/ploidetect"
    log:
        "{output_dir}/logs/mergedbed.{case}.{tumour}.log",
    shell:
        "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc})"
        "| sed 's/chr//g' > {output}" ## Cuts out any "chr" if using hg38
        " 2> {log}"
        " && ls -l {output} >> {log}"
