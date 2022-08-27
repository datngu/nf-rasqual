#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --nodes=1                
#SBATCH --job-name=vcf   
#SBATCH --mem=32G                
#SBATCH --partition=smallmem     
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
beagle4=/mnt/users/ngda/sofware/beagle.27Jan18.7e1.jar
beagle5=/mnt/users/ngda/sofware/beagle.22Jul22.46e.jar


cd /mnt/users/ngda/ngs_data/atlantic_salmon/wgs

for i in {1..29}
do
    bcftools index -t -f vcf_cp_raw/ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz
done


## phasing
mkdir vcf_phased
for i in {1..29}
do
    # estimating genotypes
    java -jar $beagle4 gtgl=vcf_cp_raw/ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz nthreads=16 niterations=20 gprobs=true out=vcf_phased/ssa${i}_tem
    # phasing
    java -jar $beagle4 gt=vcf_phased/ssa${i}_tem.vcf.gz nthreads=16 niterations=20 gprobs=true out=vcf_phased/ssa${i}
done


## add GP
for i in {1..29}
do
    /mnt/users/ngda/sofware/add_GP_vcf.py vcf_phased/ssa${i}_tem.vcf.gz vcf_phased/ssa${i}.vcf.gz vcf_phased/ssa${i}_added_GP.vcf
    bgzip -f vcf_phased/ssa${i}_added_GP.vcf
done

for i in {1..29}
do
    echo vcf_phased/ssa${i}_added_GP.vcf.gz >> file_list.txt
done

bcftools concat -n -f file_list.txt -Oz -o all_chr_added_GP.vcf.gz