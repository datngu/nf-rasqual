#!/bin/bash
in=$1
out=$2
chr=$3

# in=atac_count.tsv
# out=20.tsv
# chr=20

awk 'NR==1{print }' $in > $out
awk -v chr=$chr '{ if ($2 == chr) { print } }' $in >> $out