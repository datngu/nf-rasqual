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
params.atac_bam        = "$baseDir/data/atac_bam/*.bam"
params.atac_count      = "$baseDir/data/atac_count.tsv"
params.genotype        = "$baseDir/data/genotype.vcf.gz"
params.meta            = "$baseDir/data/meta/brain.csv"
params.outdir          = "results"

// running options
params.chrom           = 1..22 
params.peer            = 1..20 
params.genotype_PCs    = 3 
params.exp_prop        = 0.5
params.maf             = 0.05
params.fdr             = 0.1
params.atac_window     = 50000
params.eqtl_window     = 1000000

// pipeline options
params.attac_qtl         = true
params.eqtl_qtl          = true
// params.expression_qtl    = true
// params.circexplorer2     = true
// params.salmon            = true
// params.leafcutter        = true

log.info """\
================================================================
                        nf-rasqual
================================================================
    genome              : $params.genome
    cdna                : $params.cdna
    annotation          : $params.annotation
    atac_bam            : $params.atac_bam
    atac_count          : $params.atac_count
    genotype            : $params.genotype 
    meta                : $params.meta
    outdir              : $params.outdir
    chrom               : $params.chrom
    peer                : $params.peer
    maf                 : $params.maf
    fdr                 : $params.fdr
    eqtl_window         : $params.eqtl_window
    atac_window         : $params.atac_window
    genotype_PCs        : $params.genotype_PCs
    attac_qtl           : $params.attac_qtl
    eqtl_qtl            : $params.eqtl_qtl
================================================================
"""

nextflow.enable.dsl=2


workflow {

    /// channel general processing
    atac_bam_ch = channel.fromPath( params.atac_bam, checkIfExists: true )
    chrom_list_ch = channel.from(params.chrom)
    peer_list_ch = channel.from(params.peer)

    /// ATAC QTL
    //atac_bam_ch.collect().view()
    BAM_rename(params.meta, atac_bam_ch.collect())
    //BAM_rename.out.view()
    ADD_AS_vcf()
}


process BAM_rename {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'bam_dir'
    memory '8 GB'

    input:
    path meta
    path bamfiles

    output:
    path "copied_files/*{.bam,.bai}"

    script:
    """
    rename_bam.R ${meta}
    """
}


process ADD_AS_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'AS_vcf'
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    path "genotype_added_AS.vcf.gz"

    script:
    """
    atac_add_AS_vcf.sh $in_vcf genotype_added_AS.vcf.gz
    """
}