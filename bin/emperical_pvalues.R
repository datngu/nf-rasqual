#!/usr/bin/env Rscript

# Author: Dat T Nguyen <ndat@utexas.edu>
# Date: 30 Sep 2022

require(data.table)
#require(qvalue)

options(stringsAsFactors=FALSE)
syntax='Usage:
              ./this_script.R out_file rasqual_result rasqual_permute_1 rasqual_permute_2 rasqual_permute_3 ...'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax, useSource = TRUE)
  quit()
}

out_file = args[1]
rasqual_result = args[2]
in_files = args[-c(1,2)]

cat("\nout file: ", out_file, "\n")
cat("\nrasqual_result: ", rasqual_result, "\n")
x = paste(in_files, collapse = " ,")
cat("\nrasqual permutation list:\n", x, "\n")


# rasqual_res_path = "permute_10_all_chromosome_rasqual_lead_snp.txt"
get_pvalue <- function(rasqual_res_path){
  d = read.delim(rasqual_res_path, header = F)
  c = d$V11
  p = pchisq(c,1)
  return(p)
}

#rasqual_result = "/net/fs-2/scale/OrionStore/Scratch/ngda/paper1/nf-rasqual/results_rasqual/all_chromosome_rasqual_lead_snp.txt"
#in_files <- list.files(pattern = "^permute_")

df = read.delim(rasqual_result, header = F)
res = df[,c(1,2)]
colnames(res) = c("feature_id", "snp_id")
res$test_pvalue = pchisq(df$V11,1)
for(file in in_files){
  col = gsub("_all_chromosome_rasqual_lead_snp.txt", "", file)
  res[,col] = get_pvalue(file)
}