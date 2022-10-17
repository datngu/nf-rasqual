#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./RNA_filtering.R in_count in_tpm out_count exp_prop'



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

in_count = args[1]
in_tpm = args[1]
out_fn = args[3]
exp_prop = as.integer(args[4])


require(data.table)

# setwd("/mnt/SCRATCH/ngda/nf-rasqual")

# in_count = "/mnt/users/ngda/ngs_data/atlantic_salmon/brain/rna_gene_level_count_salmon.txt"
# in_tpm = "/mnt/users/ngda/ngs_data/atlantic_salmon/brain/rna_gene_level_tpm_salmon.txt"
# out_fn = "rna_gene_level_count_salmon_filtered.txt"
# exp_prop = 0.5

count = fread(in_count, header = T)
tpm = fread(in_tpm, header = T)

tpm2 = tpm[,-c(1,2)]
tpm3 = tpm2 > 0.1
pick = rowSums(tpm3)/ncol(tpm3) >= exp_prop
pick_gene = tpm[[1]][pick]

count2 = count[count[[1]] %in% pick_gene,]

fwrite(count2, file = out_fn, sep = "\t")