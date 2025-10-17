#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=Brain   
#SBATCH --mem=4G                
#SBATCH --mail-user=n.dat@outlook.com
#SBATCH --mail-type=ALL


module load Nextflow/21.03
module load singularity/rpm


export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA3OTAxfS4xNGY5NTFmOWNiZmEwNjZhOGFkYzliZTg3MDc4YWI4ZTRiYTk4ZmI5"


export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity


# download reference genome and annotation
mkdir -p data/ref
curl -Lk https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz \
  -o data/ref/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz
curl -Lk https://ftp.ensembl.org/pub/release-106/gtf/salmo_salar/Salmo_salar.Ssal_v3.1.106.gtf.gz \
  -o data/ref/Salmo_salar.Ssal_v3.1.106.gtf.gz
gunzip -k data/ref/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz
gunzip -k data/ref/Salmo_salar.Ssal_v3.1.106.gtf.gz

# download genotype reference
curl -Lk https://salmobase.org/datafiles/datasets/nf-rasqual-data/ALL_chrome_phased_filtered_HWE_1e6_biSNPs_MAF_0.01.vcf.gz \
  -o data/ref/ALL_chrome_phased_filtered_HWE_1e6_biSNPs_MAF_0.01.vcf.gz

# download ATAC

bash download_atac.sh Brain

# ## download RNA data
bash download_rna.sh Brain

## results directory
mkdir -p results/Brain


## run nextflow for Brain ATAC + RNA
nextflow run main.nf -profile orion,singularity \
  --meta "$PWD/data/meta/Brain.csv" \
  --genome "$PWD/data/ref/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa" \
  --annotation "$PWD/data/ref/Salmo_salar.Ssal_v3.1.106.gtf" \
  --atac_bam "$PWD/data/Brain/atac_bam/*.bam" \
  --atac_count "$PWD/data/Brain/atac_consensus_peak_featureCounts.txt" \
  --rna_bam "$PWD/data/Brain/rna_bam/*.bam" \
  --rna_count "$PWD/data/Brain/rna_gene_level_count_salmon.txt" \
  --genotype "$PWD/data/genotype.vcf.gz" \
  --outdir "$PWD/results/Brain" \
  --atac_window 10000 \
  --external_ld true \
  --ld_genotype "$PWD/data/ref/ALL_chrome_phased_filtered_HWE_1e6_biSNPs_MAF_0.01.vcf.gz" \
  --atac_qtl true \
  --eqtl_qtl true