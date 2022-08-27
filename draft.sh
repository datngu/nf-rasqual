# aqua faang dir
cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC


module load BCFtools/1.10.2-GCC-8.3.0


###############################
# working dir
cd /mnt/users/ngda/proj/paper1

atac_bam=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/results/bwa/mergedLibrary
atac_count=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt
genotype=/mnt/users/ngda/ngs_data/atlantic_salmon/wgs/all_chr_added_GP.vcf.gz

nextflow run main.nf -resume --atac_bam $atac_bam --atac_count $atac_count --genotype $genotype --outdir $outdir