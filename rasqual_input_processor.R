

#setwd("/Users/datn/github/nf-rasqual/data")
require(rasqualTools)
require(data.table)


# author Natsuhiko Kumasaka
randomize <- function(x,g=NULL){
  if(is.null(g)){
    n=ncol(x);
    t(apply(x,1,function(xx){xx[order(runif(n))]}))
  }else{
    for(i in unique(g)){
      x[,g==i]=randomize(x[,g==i,drop=F])
    }
    x
  }
}

rasqualMakeCovariates <- function(counts, size_factors) {
  
  #Map parameters to Natsuhiko's variables
  Y = counts
  K = size_factors
  n=ncol(Y)
  
  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9
  
  # Singular value decomposition
  fpkm.svd   = svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
  fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))
  
  # Covariate selection
  sf=log(apply(Y,2,sum))
  covs=fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
  if(cor(sf,covs[,1])^2<0.9){covs=cbind(sf, covs)}
  
  # Write covariates
  return(covs)
}


meta_fn = "meta/brain.csv"
count_fn = "atac_count.tsv"
geno_fn = "genotype.vcf.gz"

meta = fread(meta_fn, sep = ",")
meta = as.data.frame(meta)
count = fread(count_fn, skip = 1, sep = "\t")
count = as.data.frame(count)

genotype = fread(geno_fn, skip = "CHROM", sep = "\t")
genotype = as.data.frame(genotype)

# ordering count data by genotype data
geno_id = colnames(genotype)[-c(1:9)]
meta = meta[meta$genotype_id %in% geno_id,]
od = match(geno_id, meta$genotype_id)
meta = meta[od,]

atac_peaks = paste(count$Geneid, count$Chr, count$Start, count$End, count$Strand, count$Length, sep = ":")
count2 = count[,-c(1:6)]
row.names(count2) = atac_peaks
count2 = count2[ ,meta$atac_count_id]
colnames(count2) = meta$genotype_id

## save count maxtrix
saveRasqualMatrices(list( atac = count2), ".", file_suffix = "exp")

## size factor
size_factors = rasqualCalculateSampleOffsets(count2, gc_correct = FALSE)
saveRasqualMatrices(list(atac = size_factors), ".", file_suffix = "size_factors")

## covariates
covs = rasqualMakeCovariates(count2,size_factors)
covs = cbind(meta[,-c(1:6)], covs)
saveRasqualMatrices(list(atac = covs), ".", file_suffix = "covs")

## counting snps overlapping atac peaks
peak_info = data.frame(gene_id = atac_peaks)
peak_info$chr = as.character(count$Chr)
peak_info$strand = 1
peak_info$strand[which(count$Strand == "-")] = -1
peak_info$strand = as.integer(peak_info$strand)
peak_info$exon_starts = as.character(count$Start)
peak_info$exon_ends = as.character(count$End)

snp_info = genotype[,c(1:3)]
colnames(snp_info) = c("chr", "pos", "snp_id")
snp_info$chr = gsub("ssa0", "", snp_info$chr, fixed = TRUE)
snp_info$chr = gsub("ssa", "", snp_info$chr, fixed = TRUE)

snp_counts = countSnpsOverlapingExons(peak_info, snp_info, cis_window = 5e5)






