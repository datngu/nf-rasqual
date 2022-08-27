# aqua faang dir
cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC




###############################
# working dir
module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load nextflow/20.10.0
module load singularity/rpm


cd /mnt/users/ngda/proj/paper1/nf-rasqual
git pull

atac_bam=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/results/bwa/mergedLibrary/*.bam*
atac_count=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt
genotype=/mnt/users/ngda/ngs_data/atlantic_salmon/wgs/all_chr_added_GP.vcf.gz

outdir=results

nextflow run main.nf -resume --atac_bam $atac_bam --atac_count $atac_count --genotype $genotype --outdir $outdir