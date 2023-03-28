#! /bin/bash

set -eu 
set -o pipefail

echo -n "please type the rank of the sample : A=NNN...CCC...  B=NCNC.... "
read augument

#to run prep_CIRIquant

prep_CIRIquant -i sample.list \
               --lib library_info.csv \
               --bsj circRNA_bsj.csv \
               --circ circRNA_info.csv \
               --ratio circRNA_ratio.csv

prepDE.py -i sample_gene.list

# filter read count (countain  2 reads in each group of samples )
# Control, Control,..., tumor, tumor, ...

[[ $augument == "A" ]] &&
cat circRNA_bsj.csv| \
awk 'BEGIN{OFS=FS=","}NR==1{print $0}NR>1{a=0;b=0;for(i=2;i<=(NF+1)/2;i++){if($i>=2){a++}};for(j=(NF+1)/2+1;j<=NF;j++){if($j>=2){b++}};if(a>=(NF-1)/2 || b>=(NF-1)/2){print $0}}' >circRNA_bsj_filter.csv \

# Control, tumor, Control, tumor, ...
[[ $augument == "B" ]] &&
cat circRNA_bsj.csv| \
awk 'BEGIN{OFS=FS=","}NR==1{print $0}NR>1{a=0;b=0;for(i=2;i<=NF;i+=2){if($i>=2){a++}};for(j=3;j<=NF;j+=2){if($j>=2){b++}};if(a>=(NF-1)/2 || b>=(NF-1)/2){print $0}}' >circRNA_bsj_filter.csv


CIRI_DE_replicate --lib library_info.csv \
                  --bsj circRNA_bsj_filter.csv \
                  --gene gene_count_matrix.csv \
                  --out circRNA_de.tsv


join -t $'\t' -o auto -a 1  <(cat circRNA_de.tsv|sed -e 1d -e 's/,/\t/g'|sort) <(cat circRNA_info.csv|sed -e 1d -e 's/^"//g' -e 's/"$//g' -e 's/","/\t/g'|sort)| \
join -t $'\t' -o auto -a 1 - <(cat circRNA_bsj.csv|sed -e 1d -e 's/,/\t/g'|sort) > final_merged_result.tsv
