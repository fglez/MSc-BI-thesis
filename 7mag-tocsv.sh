#!/bin/bash

prefix=$1

scp mcontreras@148.247.230.5:/LUSTRE/usuario/mcontreras/thesisseq/cult2021/${prefix}/checkm/storage/bin_stats_ext.tsv .

cp bin_stats_ext.tsv ${prefix}checkm.csv #make a copy of the checkm report to reformat it to readable csv file

sed -i 's/[ \t]/ /g' ${prefix}checkm.csv

sed -i 's/ {/, {/g' ${prefix}checkm.csv

sed -i 's/{//g' ${prefix}checkm.csv

sed -i 's/}//g' ${prefix}checkm.csv

cut -d, -f1,12,13 ${prefix}checkm.csv > small${prefix}check.csv #make a small version with only bin name, completeness and contamination

rm ${prefix}checkm.csv #if you want to remove the big csv file remove the hashtag
