#!/bin/bash
#
# get_allele_freqs.bash
# 1. Performs variant calling on the provided normal.bam and filters the vcf for heterozygous snps
# 2. Computes beta-allele frequences at the heterozygous snp loci in the tumor.bam
#
# Usage function
function usage()
{
    echo "get_allele_freqs.bash"
    echo ""
    echo "usage:"
    echo ""
    echo "get_allele_freqs.bash normal.bam tumour.bam hg19.fa positions.txt tmp_dir/"
    echo ""
    echo "Variables:"
    echo "normal.bam - The germline bam file"
    echo "tumour.bam - The tumour bam file"
    echo "hg19.fa - hg19 genome fasta file, must be in the same directory as the fai"
    echo "positions.txt - postions of germline SNPs, two column file formatted as chr pos"
    echo "out_dir/ - directory where files can go"
}
#
# If no args given, print usage and exit
if [ "$1" == "" ]
then
    usage
    exit
fi
#
# Create directory if doesn't exist
[[ -d $5 ]] || mkdir $5
#
# First get vcf of all heterozygous germline SNPs
samtools mpileup $1 -l $4 -f $3 -v -B | bcftools call -c | grep "0/1" > $5/het_snps.vcf
#
# Get positions of heterozygous snps for mpileup
cut -f1,2 $5/het_snps.vcf > $5/het_snp_positions.txt
#
# Now get allele frequencies in somatic using samtools mpileup - bcftools doesn't give the needed information
# To-do: Add line breaks to the awk "one-liner"
samtools mpileup $2 -l $5/het_snp_positions.txt -f $3 -B | awk -v OFS="\t" 'BEGIN{a=0;t=0;c=0;g=0;wt=0}{seq=tolower($5);a=gsub("a", "a", seq);t=gsub("t", "t", seq);c=gsub("c","c",seq);g=gsub(/g/,"",seq);wt=gsub("\\.|\\,", "", seq);print $1, $2, $3, $4, a, t, c, g, wt}' | cat <(echo -e "chr\tpos\tref\tcov\ta\tt\tc\tg\twt") - | awk 'NR == 1 {sep="\t";for (i = 5; i <= 8; i++) bases[i] = $i;next};{sep="\t";max = $5;max_ind=5;for (i = 5; i <= 8; i++) if ($i > max) {max = $i; max_ind=i};printf "%s%s%s%s%s%s%s%s%s%s%s", $1, sep, $2, sep, toupper($3), sep, $9, sep, toupper(bases[max_ind]), sep, max;printf "%s", "\n"}' > $5/loh_raw.txt
