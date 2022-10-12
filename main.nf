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
params.atac_count      = "$baseDir/data/atac_consensus_peak_featureCounts.txt"
params.rna_bam         = "$baseDir/data/rna_bam/*.bam"
params.rna_count       = "$baseDir/data/rna_gene_level_count_salmon.txt"
params.genotype        = "$baseDir/data/genotype.vcf.gz"
params.meta            = "$baseDir/data/meta/brain.csv"
params.outdir          = "results"

// running options
params.chrom           = 1..29 
params.permute         = 20
params.phenotype_PCs   = 2 
params.exp_prop        = 0.5
params.maf             = 0.05
params.fdr             = 0.1
params.atac_window     = 10000
params.eqtl_window     = 500000

// pipeline options
params.atac_qtl          = true
params.eqtl_qtl          = true


log.info """\
================================================================
                        nf-rasqual
================================================================
    genome              : $params.genome
    cdna                : $params.cdna
    annotation          : $params.annotation
    atac_bam            : $params.atac_bam
    atac_count          : $params.atac_count
    rna_bam             : $params.rna_bam
    rna_count           : $params.rna_count
    genotype            : $params.genotype 
    meta                : $params.meta
    outdir              : $params.outdir
    chrom               : $params.chrom
    permute             : $params.permute
    maf                 : $params.maf
    fdr                 : $params.fdr
    eqtl_window         : $params.eqtl_window
    atac_window         : $params.atac_window
    phenotype_PCs       : $params.phenotype_PCs
    atac_qtl            : $params.atac_qtl
    eqtl_qtl            : $params.eqtl_qtl
================================================================
"""

nextflow.enable.dsl=2


workflow {

    // channel general processing
    chrom_list_ch = channel.from(params.chrom)
    //chrom_list_ch.collect().toList().view()
    // ATAC QTL
    if( params.atac_qtl ){
        atac_bam_ch = channel.fromPath( params.atac_bam, checkIfExists: true )
        ATAC_BAM_rename(params.meta, atac_bam_ch.collect())
        ATAC_ADD_AS_vcf(params.genotype, ATAC_BAM_rename.out)

        ATAC_PROCESS_covariates(params.meta, params.atac_count, params.genotype)
        ATAC_SPLIT_chromosome(chrom_list_ch, ATAC_ADD_AS_vcf.out, params.atac_count )
        ATAC_PREPROCESS_rasqual(chrom_list_ch, params.meta, ATAC_SPLIT_chromosome.out.collect(), params.genome)

        ATAC_RUN_rasqual(chrom_list_ch, ATAC_PREPROCESS_rasqual.out.collect(), ATAC_SPLIT_chromosome.out.collect(), ATAC_PROCESS_covariates.out)
        ATAC_RUN_rasqual_permutation(chrom_list_ch, ATAC_PREPROCESS_rasqual.out.collect(), ATAC_SPLIT_chromosome.out.collect(), ATAC_PROCESS_covariates.out)

        ATAC_MERGE_rasqual(chrom_list_ch.max(), ATAC_RUN_rasqual.out.collect())
        ATAC_MERGE_rasqual_permutation(chrom_list_ch.max(), ATAC_RUN_rasqual_permutation.out.collect())
        ATAC_COMPUTE_rasqual_emperical_pvalues(ATAC_MERGE_rasqual.collect(), ATAC_MERGE_rasqual_permutation.collect())
    }

    if( params.eqtl_qtl ){
        rna_bam_ch = channel.fromPath( params.rna_bam, checkIfExists: true )
        RNA_BAM_rename(params.meta, rna_bam_ch.collect())
        RNA_ADD_AS_vcf(params.genotype, RNA_BAM_rename.out)
        RNA_SPLIT_chromosome(chrom_list_ch, RNA_ADD_AS_vcf.out, params.rna_count )
    }
}

// rename BAM

process ATAC_BAM_rename {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_bam_dir', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path bamfiles

    output:
    path "copied_files/*{.bam,.bai}"

    script:
    """
    ATAC_rename_bam.R ${meta}
    """
}

process RNA_BAM_rename {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'RNA_bam_dir', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path bamfiles

    output:
    path "copied_files/*{.bam,.bai}"

    script:
    """
    RNA_rename_bam.R ${meta}
    """
}




// add allel specific inforation

process ATAC_ADD_AS_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_AS_vcf', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    tuple path("processed.vcf.gz"), path("processed.vcf.gz.tbi")

    script:
    """
    ls \$PWD/*bam > bam_list.txt
    zcat $in_vcf | sed 's/ssa0//g' | sed 's/ssa//g' | bgzip > tem.vcf.gz
    bcftools index -t tem.vcf.gz
    createASVCF_fixed_path.sh paired_end bam_list.txt tem.vcf.gz processed.vcf.gz rna
    bcftools index -t processed.vcf.gz
    rm tem.vcf.gz tem.vcf.gz.tbi
    """
}


process RNA_ADD_AS_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'RNA_AS_vcf', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    tuple path("processed.vcf.gz"), path("processed.vcf.gz.tbi")

    script:
    """
    ls \$PWD/*bam > bam_list.txt
    zcat $in_vcf | sed 's/ssa0//g' | sed 's/ssa//g' | bgzip > tem.vcf.gz
    bcftools index -t tem.vcf.gz
    createASVCF_fixed_path.sh paired_end bam_list.txt tem.vcf.gz processed.vcf.gz atac
    bcftools index -t processed.vcf.gz
    rm tem.vcf.gz tem.vcf.gz.tbi
    """
}



// PCA

process ATAC_PROCESS_covariates {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_covariates', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path atac_count
    path in_vcf

    output:
    tuple path("atac.covs_all_chrom.bin"), path("atac.covs_all_chrom.txt")

    script:
    """
    ATAC_covariates.R $meta $atac_count $in_vcf $params.phenotype_PCs
    """
}



// slipt chomosome

process ATAC_SPLIT_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_split_chrom', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path in_vcf
    path in_exp

    output:
    tuple path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), path("${chr}_count.txt")

    script:
    """
    awk 'NR==1{print }' $in_exp > ${chr}_count.txt
    awk 'NR==2{print }' $in_exp >> ${chr}_count.txt
    awk -v chr=$chr '{ if (\$2 == $chr) { print } }' $in_exp >> ${chr}_count.txt

    bcftools view processed.vcf.gz --regions $chr -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    """
}


process RNA_SPLIT_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'RNA_split_chrom', mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path in_vcf
    path in_exp

    output:
    tuple path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), path("${chr}_count.txt")

    script:
    """
    awk 'NR==1{print }' $in_exp > ${chr}_count.txt
    awk -v chr=$chr '{ if (\$2 == $chr) { print } }' $in_exp >> ${chr}_count.txt
    
    bcftools view processed.vcf.gz --regions $chr -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    """
}










// proprocessing


process ATAC_PREPROCESS_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_qtl_input', mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 8

    input:
    val chr
    path meta
    path split_chrom
    path genome

    output:
    tuple path("${chr}_atac.exp.bin"), path("${chr}_atac.exp.txt"), path("${chr}_atac.size_factors.bin"), path("${chr}_atac.size_factors.txt"), path("${chr}_snp_counts.tsv")


    script:
    """
    ATAC_rasqual_processor.R ${meta} ${chr}_count.txt ${chr}.vcf.gz $genome $params.atac_window ${task.cpus}
    ## rename files
    mv atac.exp.bin ${chr}_atac.exp.bin
    mv atac.exp.txt ${chr}_atac.exp.txt
    mv atac.size_factors.bin ${chr}_atac.size_factors.bin
    mv atac.size_factors.txt ${chr}_atac.size_factors.txt
    mv snp_counts.tsv ${chr}_snp_counts.tsv
    """
}






process ATAC_RUN_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_results_rasqual', mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path preproces_data
    path split_chrom
    path covariates

    output:
    path("${chr}_rasqual_lead_snp.txt")


    script:
    """
    rasqual.R vcf=${chr}.vcf.gz y=${chr}_atac.exp.bin k=${chr}_atac.size_factors.bin x=atac.covs_all_chrom.bin x_txt=atac.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_lead_snp.txt cpu=${task.cpus}
    """
}


process ATAC_MERGE_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_results_rasqual', mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    path atac_rasqual_results

    output:
    path("all_chromosome_rasqual_lead_snp.txt")


    script:
    """
    for chr in \$(seq 1 $max_chr)
    do
        cat \${chr}_rasqual_lead_snp.txt >> all_chromosome_rasqual_lead_snp.txt
    done
    """
}


process ATAC_RUN_rasqual_permutation {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_results_rasqual_permutaion', mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path preproces_data
    path split_chrom
    path covariates

    output:
    path("${chr}_permute_*_rasqual_lead_snp.txt")


    script:
    """
    rasqual_permute.sh ${chr}.vcf.gz ${chr}_atac.exp.bin ${chr}_atac.size_factors.bin atac.covs_all_chrom.bin atac.covs_all_chrom.txt ${chr}_snp_counts.tsv ${task.cpus} ${params.permute} ${chr}
    """
}


process ATAC_MERGE_rasqual_permutation {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_results_rasqual_permutaion', mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    path atac_rasqual_results

    output:
    path("permute_*_all_chromosome_rasqual_lead_snp.txt")


    script:
    """
    for i in \$(seq 1 $params.permute)
    do
        for chr in \$(seq 1 $max_chr)
        do
            cat \${chr}_permute_\${i}_rasqual_lead_snp.txt >> permute_\${i}_all_chromosome_rasqual_lead_snp.txt
        done
    done
    """
}



process ATAC_COMPUTE_rasqual_emperical_pvalues {
    container 'ndatth/rasqual:v0.0.0'
    publishDir 'ATAC_results_emperical_pvalues', mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path merged_results
    path permuation_merged_results

    output:
    path("rasqual_emperical_pvalues.txt")


    script:
    """
    rasqual_emperical_pvalues.R rasqual_emperical_pvalues.txt $merged_results $permuation_merged_results
    """
}