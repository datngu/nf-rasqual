#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-rasqual
========================================================================================
    RASQUAL Analysis Pipeline.
 https://github.com/datngu/nf-rasqual
 Author: Dat T Nguyen
 Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/





/*
 Define the default parameters
*/ 
params.genome          = "$baseDir/data/ref/genome.fa"
params.cdna            = "$baseDir/data/ref/cdna.fa"
params.annotation      = "$baseDir/data/ref/annotation.gtf"
params.atac_reads      = "$baseDir/data/atac_reads/*_{1,2}.fastq.gz"
params.rna_reads       = "$baseDir/data/rna_reads/*_{1,2}.fastq.gz"
params.genotype        = "$baseDir/data/genotype.vcf.gz"
params.meta            = "$baseDir/data/meta.csv"
params.outdir          = "results"

// running options
params.chrom           = 1..22 
params.peer            = 1..20 
params.genotype_PCs    = 4 
params.exp_prop        = 0.5
params.fdr             = 0.1
params.cis_window      = 500000

// pipeline options
params.attac_qtl         = true
params.expression_qtl    = true
params.circexplorer2     = true
params.salmon            = true
params.leafcutter        = true

log.info """\
================================================================
                        nf-rasqual
================================================================
    genome              : $params.genome
    cdna                : $params.cdna
    bsj                 : $params.bsj
    annotation          : $params.annotation
    reads               : $params.reads
    genotype            : $params.genotype 
    meta                : $params.meta
    outdir              : $params.outdir
    chrom               : $params.chrom
    bsj_filter          : $params.bsj_filter
    exp_prop            : $params.exp_prop
    peer                : $params.peer
    maf                 : $params.maf
    fdr                 : $params.fdr
    fastqtl_window      : $params.fastqtl_window
    genotype_PCs        : $params.genotype_PCs
    circall             : $params.circall
    ciri2               : $params.ciri2
    circexplorer2       : $params.circexplorer2
    salmon              : $params.salmon
    leafcutter          : $params.leafcutter
================================================================
"""


workflow {



}