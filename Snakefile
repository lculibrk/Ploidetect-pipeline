configfile: "./CONFIG.txt"
chromosomes=config["chromosomes"]

cases=config["bams"].keys()

rule all:
    input:
        expand("data/{id}/segmented.RDS", id = cases)
    
rule germline_cov:
    input:
            lambda wildcards: config["bams"][wildcards.id]["normal"]
    output:
            "data/{id}/normal/{chr}.bed"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "samtools depth -r{wildcards.chr} -Q 21 {input} | awk -v FS='\t' -v OFS='\t' 'NR > 1{{print $1, $2, $2+1, $3}}' | python3 scripts/make_windows.py - 100000 | bedtools sort -i stdin > {output}"

rule merge_germline:
    input:
            expand("data/{{id}}/normal/{chr}.bed", chr=chromosomes)
    output:
            "data/{id}/germline.bed"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "cat {input} | bedtools sort -i stdin > {output}"

rule makewindowfile:
    input:
            rules.merge_germline.output
    output:
            "data/{id}/windows.txt"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "cut -f 1,2,3 < {input} | bedtools sort -i stdin > {output}"
rule splitwindowfile:
    input:
            rules.makewindowfile.output
    output:
            "data/{id}/windows/{chr}.txt"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"
rule genomecovsomatic:
    input:
            lambda wildcards: config["bams"][wildcards.id]["somatic"],
            lambda wildcards: config["bams"][wildcards.id]["normal"],
            window=rules.splitwindowfile.output
    output:
            "data/{id}/tumour/{chr}.bed"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "bedtools multicov -bams {input[0]} {input[1]} -q 20 -bed {input.window}  > {output}"

rule mergesomatic:
    input:
            expand("data/{{id}}/tumour/{chr}.bed", chr=chromosomes)
    output:
            "data/{id}/tumour.bed"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "cat {input} | bedtools sort -i stdin > {output}"

rule compute_loh:
    input:
            normbam=lambda wildcards: config["bams"][wildcards.id]["normal"],
            sombam=lambda wildcards: config["bams"][wildcards.id]["somatic"]
    output:
            "data/{id}/loh_tmp/loh_raw.txt"
    params:
            genome=config["genome"]
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "bash scripts/get_allele_freqs.bash {input.normbam} {input.sombam} {params.genome} resources/SNP_array_positions.txt data/{wildcards.id}/loh_tmp/"

rule process_loh:
    input:
            loh=rules.compute_loh.output,
            window=rules.makewindowfile.output
    output:
            "data/{id}/loh.bed"
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}' {input.loh} | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}' | bedtools sort -i stdin | Rscript scripts/merge_loh.R -l STDIN -w {input.window} -o {output}"

rule getgc:
    input:
                window=rules.makewindowfile.output
    output:
                "data/{id}/gc.bed"
    params:
                genome=config["genome"]
    conda:
                "conda_configs/sequence_processing.yaml"
    shell:
                "bedtools nuc -fi {params.genome} -bed {input} | cut -f1,2,3,5 | tail -n +2 > {output}"

rule mergedbed:
    input:
                gc=rules.getgc.output,
                tumour=rules.mergesomatic.output,
                normal=rules.merge_germline.output,
		loh=rules.process_loh.output
    output:
                "data/{id}/merged.bed"
    conda:
                "conda_configs/sequence_processing.yaml"
    shell:
                "paste {input.tumour} <(cut -f4 {input.loh}) <(cut -f4 {input.gc}) > {output}"

rule preseg:
    input:
        "data/{id}/merged.bed"
    output:
        "data/{id}/segmented.RDS"
    shell:
        "Rscript scripts/prep_ploidetect2.R -i {input} -o {output}"
