import glob
import os
configfile: "./CONFIG.txt"
chromosomes=config["chromosomes"]
genome=config["genome_name"]

array_positions = glob.glob(os.path.join("resources/snp_arrays", genome, "SNP_array_positions.txt"))[0]

print(array_positions)

cases=config["bams"].keys()

rule all:
    input:
        expand("output/{id}/cna.txt", id = cases)
    
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
            "bash scripts/get_allele_freqs.bash {input.normbam} {input.sombam} {params.genome} " + array_positions +  " data/{wildcards.id}/loh_tmp/"

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

rule install_pdt:
    output:
        "pdt_installed.txt"
    conda:
        "conda_configs/r.yaml"
    shell:
        "Rscript -e \"devtools::install_github('lculibrk/Ploidetect')\"; touch {output}"

rule preseg:
    input:
        "data/{id}/merged.bed"
    output:
        "data/{id}/segmented.RDS"
    conda:
        "conda_configs/r.yaml"
    shell:
        "Rscript scripts/prep_ploidetect2.R -i {input} -o {output}"

rule ploidetect:
    input:
        rules.preseg.output,
	rules.install_pdt.output
    output:
        "output/{id}/plots.pdf",
        "output/{id}/models.txt",
        "output/{id}/meta.RDS"
    conda:
        "conda_configs/r.yaml"
    resources: cpus=24, mem_mb=189600
    shell:
        "Rscript scripts/run_ploidetect2.R -i {input[0]} -o {output[1]} -p {output[0]} -r {output[2]}"

rule ploidetect_copynumber:
    input:
        "output/{id}/models.txt",
        "data/{id}/segmented.RDS",
        "output/{id}/plots.pdf"
    output:
        "output/{id}/cna.txt",
        "output/{id}/cna_plots.pdf"
    conda:
        "conda_configs/r.yaml"
    log:
        "output/{id}/cna_log.txt"
    resources: cpus=24, mem_mb=189600
    shell:
        "Rscript scripts/ploidetect_copynumber.R -i {input[1]} -m {input[0]} -p {output[1]} -o {output[0]} &> {log}"
												    