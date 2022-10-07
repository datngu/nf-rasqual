#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=cpfiles   
#SBATCH --mem=4G                
#SBATCH --partition=smallmem     
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL




home_dir=/mnt/users/ngda/ngs_data/atlantic_salmon/brain

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/Brain/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done


