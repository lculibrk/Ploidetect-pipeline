BEGIN{
    sep="\t"
    OFS="\t"
    a=0;t=0;c=0;g=0;wt=0;
    bases[5]="a"
    bases[6]="t"
    bases[7]="c"
    bases[8]="g"
    #print "chr", "pos", "ref", "cov", "a", "t", "c", "g", "wt"
}
{
    seq=tolower($5);
    counts[5]=gsub("a", "a", seq);
    counts[6]=gsub("t", "t", seq);
    counts[7]=gsub("c","c",seq);
    counts[8]=gsub(/g/,"",seq);
    wt=gsub("\\.|\\,", "", seq)
    max=counts[5]
    max_ind = 5;
    for (i = 5; i <= 8; i++){
        if (counts[i] > max) {
            max = counts[i]; max_ind=i;
        }
    }
    printf "%s%s%s%s%s%s%s%s%s%s%s", $1, sep, $2, sep, toupper($3), sep, wt, sep, toupper(bases[max_ind]), sep, max;
    printf "%s", "\n"
}
