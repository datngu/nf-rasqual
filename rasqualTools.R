#devtools::install_github("kauralasoo/rasqual/rasqualTools")
setwd("/Users/datn/github/nf-rasqual")
require(rasqualTools)
head(counts_matrix)
saveRasqualMatrices(list(cellTypeA = counts_matrix), ".", file_suffix = "expression")

size_factors = rasqualCalculateSampleOffsets(counts_matrix, gc_correct = FALSE)
saveRasqualMatrices(list(cellTypeA = size_factors), ".", file_suffix = "size_factors")

covs = rasqualMakeCovariates(counts_matrix,size_factors)

fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9

fpkm_norm = log2(fpkm+0.1)
# Singular value decomposition
x = (fpkm_norm-apply(fpkm_norm,1,mean))/apply(fpkm_norm,1,sd)
fpkm.svd   = svd(x)

fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))