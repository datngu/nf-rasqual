#devtools::install_github("kauralasoo/rasqual/rasqualTools")
setwd("/Users/datn/github/nf-rasqual")
require(rasqualTools)
head(counts_matrix)
saveRasqualMatrices(list(cellTypeA = counts_matrix), ".", file_suffix = "expression")

size_factors = rasqualCalculateSampleOffsets(counts_matrix, gc_correct = FALSE)
saveRasqualMatrices(list(cellTypeA = size_factors), ".", file_suffix = "size_factors")
