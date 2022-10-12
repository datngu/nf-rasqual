setwd("/net/fs-2/scale/OrionStore/Scratch/ngda/nf-rasqual/results_rasqual_permutaion")
rasqual_result = "/net/fs-2/scale/OrionStore/Scratch/ngda/nf-rasqual/results_rasqual/all_chromosome_rasqual_lead_snp.txt"

require(data.table)
require(qvalue)


get_pvalue <- function(rasqual_res_path){
  d = fread(rasqual_res_path, header = F, fill = T)
  c = d$V11
  p = pchisq(c,1)
  return(p)
}

get_chi_square <- function(rasqual_res_path){
  d = fread(rasqual_res_path, header = F, fill = T)
  c = d$V11
  return(c)
}

df = fread(rasqual_result, header = F, fill = T)
res = df[,c(1,2)]
colnames(res) = c("feature_id", "snp_id")
res$test_pvalue = pchisq(df$V11,1)
#res$chi_square = df$V11

for(file in in_files){
  col = gsub("_all_chromosome_rasqual_lead_snp.txt", "", file)
  res[,col] = get_pvalue(file)
  #res[,col] = get_chi_square(file)
}
#res2 = res
res = res[!res$snp_id == "SKIPPED", ]


x1 = as.numeric(res[[3]])
x2 = as.matrix(res[,4:33])
emp = empPvals(x1, x2, pool = FALSE)
emp_fdr = qvalue(emp)
q = emp_fdr$qvalues
sum(q < 0.1)














################
df = fread(rasqual_result, header = F, fill = T)
res_cq = df[,c(1,2)]
colnames(res_cq) = c("feature_id", "snp_id")
res_cq$chi_square = df$V11

for(file in in_files){
  col = gsub("_all_chromosome_rasqual_lead_snp.txt", "", file)
  res_cq[,col] = get_chi_square(file)
}
#res2 = res
res_cq = res_cq[!res_cq$snp_id == "SKIPPED", ]


x1 = as.numeric(res_cq[[3]])
x2 = as.matrix(res_cq[,4:33])
emp = empPvals(x1, x2, pool = FALSE)
emp = empPvals(x1, x2, pool = T)
hist(emp,150)
emp_fdr = qvalue(emp)
q = emp_fdr$qvalues
sum(q < 0.05)















