

setwd("/Users/datn/github/nf-rasqual/data")
require(rasqualTools)
require(data.table)



count_fn = "attac_count.tsv"
geno_fn = "genotype.vcf.gz"

count = fread(count_fn, skip = 1)
