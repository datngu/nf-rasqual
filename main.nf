#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-rasqual
========================================================================================
                RASQUAL Analysis Pipeline with nextflow.
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
params.atac_qtl          = true
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
    atac_qtl            : $params.atac_qtl
    eqtl_qtl            : $params.eqtl_qtl
================================================================
"""

nextflow.enable.dsl=2


workflow {

    /// channel general processing
    atac_bam_ch = channel.fromPath( params.atac_bam, checkIfExists: true )
    chrom_list_ch = channel.from(params.chrom)
    peer_list_ch = channel.from(params.peer)

    // spliting vcf and counting expression
    INDEX_vcf(params.genotype)
    SPLITING_chromosome(chrom_list_ch, INDEX_vcf.out, params.atac_count )

    /// ATAC QTL
    //atac_bam_ch.collect().view()
    BAM_rename(params.meta, atac_bam_ch.collect())
    ADD_AS_vcf(params.genotype, BAM_rename.out)
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
    publishDir 'atac_AS_vcf'
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    path "genotype_added_AS.vcf.gz"

    script:
    """
    ls \$PWD/*bam > bam_list.txt
    bcftools index -tf $in_vcf
    createASVCF_fixed_path.sh paired_end bam_list.txt $in_vcf genotype_added_AS.vcf.gz atac
    """
}


process ADD_AS_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'atac_AS_vcf'
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    path "genotype_added_AS.vcf.gz"

    script:
    """
    ls \$PWD/*bam > bam_list.txt
    bcftools index -tf $in_vcf
    createASVCF_fixed_path.sh paired_end bam_list.txt $in_vcf genotype_added_AS.vcf.gz atac
    """
}



process INDEX_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'atac_AS_vcf'
    memory '8 GB'

    input:
    path in_vcf

    output:
    tuple path("processed.vcf.gz"), path("processed.vcf.gz.tbi")

    script:
    """
    zcat $in_vcf | sed 's/ssa0//g' | sed 's/ssa//g' | bgzip > processed.vcf.gz
    bctools index -t processed.vcf.gz
    """
}



process SPLITING_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'split_chrom'
    memory '8 GB'

    input:
    val chr
    path in_vcf
    path in_atac_exp

    output:
    tuple val("${chr}"), path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), path("${chr}.atac_count.txt")

    script:
    """
    atac_exp_filter.py $in_atac_exp ${chr}.atac_count.txt $chr
    bcftools view processed.vcf.gz --regions $chr -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    """

}

process PREPROCESSING_atac_qtl {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'atac_qtl_input'
    memory '8 GB'

    input:
    path meta
    path splited_chrom

    output:
    tuple val("${chr}"), path("${chr}_atac.covs.bin"), path("${chr}_atac.covs.txt"), path("${chr}_atac.exp.bin"), path("${chr}_atac.exp.txt"), path("${chr}_atac.size_factors.bin"), path("${chr}_atac.size_factors.txt"), path("${chr}_snp_counts.tsv")


    script:
    """
    atac_rasqual_processor.R ${meta} ${chr}.atac_count.txt ${chr}.vcf.gz
    ## rename files
    mv atac.covs.bin ${chr}_atac.covs.bin
    mv atac.covs.txt ${chr}_atac.covs.txt
    mv atac.exp.bin ${chr}_atac.exp.bin
    mv atac.exp.txt ${chr}_atac.exp.txt
    mv atac.size_factors.bin ${chr}_atac.size_factors.bin
    mv atac.size_factors.txt ${chr}_atac.size_factors.txt
    mv snp_counts.tsv ${chr}_snp_counts.tsv
    """
}