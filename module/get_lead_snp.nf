#!/usr/bin/env nextflow

// include { foo } from './module/get_lead_snp'

nextflow.enable.dsl=2



// process to obtain lead snps

process ATAC_MERGE_rasqual_normial {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.trace_dir}/ATAC_MERGE_rasqual_normial", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path all_inputs

    output:
    path("obs_rasqual_normial_pval.txt")


    script:
    """
    cat *_formated_EigenMT.txt >> obs_rasqual_normial_pval.txt
    """
}


process ATAC_MERGE_rasqual_normial_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.trace_dir}/ATAC_MERGE_rasqual_normial_permute", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path all_inputs

    output:
    path("nul_rasqual_normial_pval.txt")


    script:
    """
    cat *_formated_EigenMT.txt >> nul_rasqual_normial_pval.txt
    """
}


process ATAC_GET_lead_SNP {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_GET_lead_SNP", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path eigenMT_merged
    path normial_merged

    output:
    path("obs_rasqual_normial_pval_lead_snp.txt")


    script:
    """
    extract_lead_SNPs_nomial_pval.py --normial $normial_merged --eigenMT $eigenMT_merged > obs_rasqual_normial_pval_lead_snp.txt
    """
}


process ATAC_GET_lead_SNP_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_GET_lead_SNP_permute", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path eigenMT_merged
    path normial_merged

    output:
    path("nul_rasqual_normial_pval_lead_snp.txt")


    script:
    """
    extract_lead_SNPs_nomial_pval.py --normial $normial_merged --eigenMT $eigenMT_merged > nul_rasqual_normial_pval_lead_snp.txt
    """
}



// process to obtain lead snps

process RNA_MERGE_rasqual_normial {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.trace_dir}/RNA_MERGE_rasqual_normial", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path all_inputs

    output:
    path("obs_rasqual_normial_pval.txt")


    script:
    """
    cat *_formated_EigenMT.txt >> obs_rasqual_normial_pval.txt
    """
}


process RNA_MERGE_rasqual_normial_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.trace_dir}/RNA_MERGE_rasqual_normial_permute", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path all_inputs

    output:
    path("nul_rasqual_normial_pval.txt")


    script:
    """
    cat *_formated_EigenMT.txt >> nul_rasqual_normial_pval.txt
    """
}


process RNA_GET_lead_SNP {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_GET_lead_SNP", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path eigenMT_merged
    path normial_merged

    output:
    path("obs_rasqual_normial_pval_lead_snp.txt")


    script:
    """
    extract_lead_SNPs_nomial_pval.py --normial $normial_merged --eigenMT $eigenMT_merged > obs_rasqual_normial_pval_lead_snp.txt
    """
}


process RNA_GET_lead_SNP_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_GET_lead_SNP_permute", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path eigenMT_merged
    path normial_merged

    output:
    path("nul_rasqual_normial_pval_lead_snp.txt")


    script:
    """
    extract_lead_SNPs_nomial_pval.py --normial $normial_merged --eigenMT $eigenMT_merged > nul_rasqual_normial_pval_lead_snp.txt
    """
}
