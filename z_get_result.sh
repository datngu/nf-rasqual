cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/results

rm -rf result_collect_permutation
mkdir -p result_collect_permutation


get_tis (){
    tis=$1
    mkdir result_collect_permutation/$tis
    cp $tis/ATAC_results_emperical_pvalues/rasqual_emperical_pvalues.txt result_collect_permutation/$tis/atac_rasqual_emperical_pvalues.txt
    cp $tis/RNA_results_emperical_pvalues/rasqual_emperical_pvalues.txt result_collect_permutation/$tis/rna_rasqual_emperical_pvalues.txt
    #cp $tis/deltaSVM_prediction/all_snp_deltaSVM.txt result_collect_permutation/$tis/all_snp_deltaSVM.txt
    cp $tis/ATAC_results_rasqual/all_chromosome_rasqual_lead_snp.txt result_collect_permutation/$tis/atac_all_chromosome_rasqual_lead_snp.txt
    cp $tis/RNA_results_rasqual/all_chromosome_rasqual_lead_snp.txt result_collect_permutation/$tis/rna_all_chromosome_rasqual_lead_snp.txt
    #cp $tis/ATAC_filtering_expression/atac_consensus_peak_featureCounts_filtered.txt result_collect_permutation/$tis/atac_consensus_peak_featureCounts_filtered.txt
}


get_tis Brain
get_tis Liver
get_tis Gonad
get_tis Muscle
get_tis Gill

get -r /mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/results/result_collect_permutation /Users/datn/DATA_ANALYSES/aqua_qtl_permutation