#!/bin/bash

in_vcf=$1
out_vcf=$2

ls $PWD/*bam > bam_list.txt

bcftools index -tf $in_vcf

createASVCF_fixed_path.sh paired_end bam_list.txt $in_vcf $out_vcf atac
