#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./ATAC_rasqual_processor.R meta_csv feature_count_txt genotype_vcf genome.fa cis_window phenotype_PCs n_core'

# defaut outputs:
# atac.covs.bin - atac.covs.txt
# atac.exp.bin - atac.exp.txt
# atac.size_factors.bin - atac.size_factors.txt
# snp_counts.tsv



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 5 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

cpu = NA
meta_fn = args[1]
count_fn = args[2]
geno_fn = args[3]
genome_fn = args[4]
cis_window = as.integer(args[5])
phenotype_PCs = as.integer(args[6])
if(length(args) >= 7) cpu = as.integer(args[7])
if (is.na(cpu)) cpu=1
cpu = as.integer(cpu)



#setwd("/Users/datn/github/nf-rasqual/data")

#############

# input
#  cis_window = 10000
#  phenotype_PCs = 2
#  meta_fn = "meta/brain.csv"
#  count_fn = "atac_consensus_peak_featureCounts.txt"
#  geno_fn = "genotype.vcf.gz"
#  genome_fn = "/Users/datn/GENOMES/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna.toplevel.fa"

#############
# output:
# atac.covs.bin - atac.covs.txt
# atac.exp.bin - atac.exp.txt
# atac.size_factors.bin - atac.size_factors.txt
# snp_counts.tsv
#############


require(rasqualTools)
require(GenomicFeatures)
require(Biostrings)
require(data.table)
require(foreach)
require(doParallel)
registerDoParallel(cores=cpu)


set.seed(2022)


PCA_Covariates <- function(counts, size_factors, n_PCs = 2) {
  # author Natsuhiko Kumasaka
  #Map parameters to Natsuhiko's variables
  Y = counts
  K = size_factors
  n=ncol(Y)
  
  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9
  
  ## author Dat T Nguyen
  p = prcomp(fpkm)
  pca = p$rotation[,1:n_PCs]
  # Write covariates
  return(pca)
}






get_GC <- function(genome, feature_info){
  #for(i in 1:nrow(feature_info)){
  seqs = foreach(i = 1:nrow(feature_info),.combine = c) %dopar% {
    row = feature_info[i,]
    chr = row$chr
    s = as.integer(strsplit(row$exon_starts, ","))
    e = as.integer(strsplit(row$exon_ends, ","))
    x = substring(genome[chr], s, e)
  }
  seqs = DNAStringSet(seqs)
  alf <- Biostrings::alphabetFrequency(seqs, as.prob=TRUE)
  gc = rowSums(alf[,c("G", "C"), drop=FALSE])
  return(gc)
}


meta = fread(meta_fn, sep = ",")
meta = as.data.frame(meta)
count = fread(count_fn, skip = 1, sep = "\t")
count = as.data.frame(count)
genotype = fread(geno_fn, skip = "CHROM", sep = "\t")
genotype = as.data.frame(genotype)

# process genome for GC counting
genome = readDNAStringSet(genome_fn)
tem = strsplit(names(genome), " ")
names(genome) = do.call("rbind", tem)[,1]


# ordering count data by genotype data
geno_id = colnames(genotype)[-c(1:9)]
meta = meta[meta$genotype_id %in% geno_id,]
od = match(geno_id, meta$genotype_id)
meta = meta[od,]
rownames(meta) = meta$genotype_id

atac_peaks = paste(count$Geneid, count$Chr, count$Start, count$End, count$Strand, count$Length, sep = ":")


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
#snp_info$chr = gsub("ssa0", "", snp_info$chr, fixed = TRUE)
#snp_info$chr = gsub("ssa", "", snp_info$chr, fixed = TRUE)

snp_counts = countSnpsOverlapingExons(peak_info, snp_info, cis_window = cis_window)
fwrite(snp_counts, file = "snp_counts.tsv", sep = "\t")



## save count maxtrix
count2 = count[,-c(1:6)]
row.names(count2) = atac_peaks
count2 = count2[ ,meta$atac_count_id]
colnames(count2) = meta$genotype_id
saveRasqualMatrices(list( atac = count2), ".", file_suffix = "exp")

## size factor
# GC counting
gc = get_GC(genome, peak_info)
gc_percentage = gc*100
GC = data.frame(gene_id = row.names(count2), percentage_gc_content = gc_percentage)
# comput size factors
size_factors = rasqualCalculateSampleOffsets(count2, GC)
saveRasqualMatrices(list(atac = size_factors), ".", file_suffix = "size_factors")

## covariates
covs = PCA_Covariates(count2, size_factors, phenotype_PCs)
covs = cbind(meta[,-c(1:6)], covs)
saveRasqualMatrices(list(atac = covs), ".", file_suffix = "covs")


