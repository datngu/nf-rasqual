#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./ATAC_filtering.R feature_count_txt out_fn chrom_list'




args = commandArgs(trailingOnly = TRUE)

if(length(args) < 2 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

require(data.table)

count_fn = args[1]
out_fn = args[2]
chrom_list = args[-c(1,2)]

df = data.frame(chrom_list = chrom_list)

fwrite(df, file = out_fn, sep = "\t")