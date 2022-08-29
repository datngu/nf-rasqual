#!/bin/bash

in_vcf=$1
out_vcf=$2

ls $PWD/*bam > bam_list.txt

bcftools sort $in_vcf -Oz -o sorted.vcf.gz
bcftools index -tf sorted.vcf.gz

createASVCF.sh paired_end bam_list.txt sorted.vcf.gz $out_vcf atac
