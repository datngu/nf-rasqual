#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

meta_path = out_file = args[1]
# meta_path = "data/meta/brain.csv"
meta = read.csv(meta_path, header = TRUE)

for( i in 1: nrow(meta)){
    # bam
    new_bam = paste0(meta$genotype_id[i], ".bam")
    old_bam = meta$attac_bam_id[i]
    cmd = paste0("cp ", old_bam, " ", new_bam)
    print(cmd)
    #system(cmd)
    # bai
    new_bai = paste0(meta$genotype_id[i], ".bam.bai")
    old_bai = paste0(meta$attac_bam_id[i], ".bai")
    cmd2 = paste0("cp ", old_bai, " ", new_bai)
    print(cmd2)
    #system(cmd2)
}


