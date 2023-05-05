cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/results

rm -rf result_collect_permutation
mkdir -p result_collect_permutation


get_tis (){
    tis=$1
    mkdir result_collect_permutation/$tis
    cp $tis/ATAC_results_emperical_pvalues/rasqual_emperical_pvalues.txt result_collect_permutation/$tis/atac_rasqual_emperical_pvalues.txt
    cp $tis/RNA_results_emperical_pvalues/rasqual_emperical_pvalues.txt result_collect_permutation/$tis/rna_rasqual_emperical_pvalues.txt
    cp $tis/deltaSVM_prediction/all_snp_deltaSVM.txt result_collect_permutation/$tis/all_snp_deltaSVM.txt
}


get_tis Brain
get_tis Liver

get -r /mnt/ScratchProjects/Aqua-Faang/dat_projects/aqua_qtl/results/result_collect_permutation /Users/datn/DATA_ANALYSES/aqua_qtl_permutation