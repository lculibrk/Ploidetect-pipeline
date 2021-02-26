import glob
import os
configfile: "./CONFIG.txt"
chromosomes=config["chromosomes"]
genome=config["genome_name"]

output_dir = config["output_dir"]
temp_dir = config["temp_dir"] if config["temp_dir"] else output_dir + "/temp"

rule all:
    input:
        expand("{output_dir}/cna.txt", output_dir=output_dir)

rule germline_cov:
    input:
            config["bams"]["normal"]
    output:
            temp("{temp_dir}/normal/{chr}.bed")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "samtools depth -r{wildcards.chr} -Q 21 {input} | awk -v FS='\t' -v OFS='\t' 'NR > 1{{print $1, $2, $2+1, $3}}' | python3 scripts/make_windows.py - 100000 | bedtools sort -i stdin > {output}"

rule merge_germline:
    input:
            expand("{temp_dir}/normal/{chr}.bed", chr=chromosomes, temp_dir=temp_dir)
    output:
            temp("{temp_dir}/germline.bed")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "cat {input} | bedtools sort -i stdin > {output}"

rule makewindowfile:
    input:
            rules.merge_germline.output
    output:
            temp("{temp_dir}/windows.txt")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "cut -f 1,2,3 < {input} | bedtools sort -i stdin > {output}"
rule splitwindowfile:
    input:
            rules.makewindowfile.output
    output:
            temp("{temp_dir}/windows/{chr}.txt")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "awk -v FS='\t' -v OFS='\t' '$1 == \"{wildcards.chr}\"{{print $0}}' {input} > {output}"
rule genomecovsomatic:
    input:
            sombam=config["bams"]["somatic"],
            nombam=config["bams"]["normal"],
            window=rules.splitwindowfile.output
    output:
            temp("{temp_dir}/tumour/{chr}.bed")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "bedtools multicov -bams {input.sombam} {input.nombam} -q 20 -bed {input.window}  > {output}"

rule mergesomatic:
    input:
            expand("{temp_dir}/tumour/{chr}.bed", chr=chromosomes, temp_dir=temp_dir)
    output:
            temp("{temp_dir}/tumour.bed")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "cat {input} | bedtools sort -i stdin > {output}"

rule compute_loh:
    input:
            normbam = config["bams"]["normal"],
            sombam = config["bams"]["somatic"]
    output:
            temp("{temp_dir}/loh_tmp/loh_raw.txt")
    params:
            genome=config["genome"],
            array_positions=config["array_positions"]
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "bash scripts/get_allele_freqs.bash {input.normbam} {input.sombam}"
            " {params.genome} {params.array_positions} {temp_dir}/loh_tmp/"

rule process_loh:
    input:
            loh=rules.compute_loh.output,
            window=rules.makewindowfile.output
    output:
            temp("{temp_dir}/loh.bed")
    conda:
            "conda_configs/sequence_processing.yaml"
    shell:
            "awk -v FS='\t' -v OFS='\t' '($4 != 0 && $6 != 0){{ print $1, $2, $2+1, $4, $6 }}' {input.loh} | awk -v FS='\t' -v OFS='\t' '{{print $1, $2, $3, ($4 / ($4 + $5)) }}' | bedtools sort -i stdin | Rscript scripts/merge_loh.R -l STDIN -w {input.window} -o {output}"

rule getgc:
    input:
                window=rules.makewindowfile.output
    output:
                temp("{temp_dir}/gc.bed")
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
                temp("{temp_dir}/merged.bed")
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
        "R REMOTES_NO_ERRORS_FROM_WARNINGS=\"true\"; Rscript -e \"devtools::install_github('lculibrk/Ploidetect')\"; touch {output}"

rule preseg:
    input:
        expand("{temp_dir}/merged.bed", temp_dir=temp_dir)
    output:
        temp(expand("{temp_dir}/segmented.RDS", temp_dir=temp_dir))
    conda:
        "conda_configs/r.yaml"
    shell:
        "Rscript scripts/prep_ploidetect2.R -i {input} -o {output}"

rule ploidetect:
    input:
        rules.preseg.output,
	rules.install_pdt.output
    output:
        "{output_dir}/plots.pdf",
        "{output_dir}/models.txt",
        "{output_dir}/meta.RDS"
    conda:
        "conda_configs/r.yaml"
    resources: cpus=24, mem_mb=189600
    shell:
        "Rscript scripts/run_ploidetect2.R -i {input[0]} -o {output[1]} -p {output[0]} -r {output[2]}"

rule ploidetect_copynumber:
    input:
        "{output_dir}/models.txt",
        expand("{temp_dir}/segmented.RDS", temp_dir=temp_dir),
        "{output_dir}/plots.pdf"
    output:
        "{output_dir}/cna.txt",
        "{output_dir}/cna_plots.pdf"
    conda:
        "conda_configs/r.yaml"
    log:
        "{output_dir}/cna_log.txt"
    resources: cpus=24, mem_mb=189600
    shell:
        "Rscript scripts/ploidetect_copynumber.R -i {input[1]} -m {input[0]} -p {output[1]} -o {output[0]} &> {log}"
