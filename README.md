# nf-rasqual Pipeline

## Overview

This pipeline processes and analyzes genomic, transcriptomic, and ATAC-seq data for QTL (Quantitative Trait Loci) analysis. It uses genotype and phenotype data to identify associations, focusing on eQTLs (expression QTLs) and ATAC-QTLs. The pipeline is designed for flexibility, allowing users to configure various parameters, including the use of external linkage disequilibrium (LD) data for multiple testing correction with EigenMT.



## 2. Installation and Dependencies

### Dependencies

* **Nextflow** (v21.10.6 or higher)  
* **Docker** or **Singularity** for containerized execution  

Since cscQTL is implemented with Nextflow (DSL2), you will need Nextflow to run it. Installation instructions are available on the Nextflow homepage: [https://www.nextflow.io/](https://www.nextflow.io/).  

All tool dependencies are included in the corresponding container (next section), so running with Docker or Singularity is sufficient.  

By default, nf-rasqual can run either with Docker on a single machine (**-profile local**) or with Singularity on an HPC using SLURM (**-profile hpc**). You may need to adjust **nextflow.config** to match your system configuration.  

For further customization, consult the Nextflow documentation: [https://www.nextflow.io/docs/latest/config.html](https://www.nextflow.io/docs/latest/config.html).


### Container

A pre-built Docker container includes all required software (RASQUAL, R scripts, Python scripts, bcftools, etc.) that will be automaticly downloaded  by **Nextflow**

```
ndatth/rasqual:v0.0.0
```

Alternatively, if interest, users can build the container themselves as following the instructions and using the source code and libraries provided in the `container` directory of the repository: [https://github.com/datngu/nf-rasqual/tree/main/container](https://github.com/datngu/nf-rasqual/tree/main/container).


### Hardware Requirements

The most resource-intensive step of the pipeline requires **64 GB of RAM** and **16 CPU cores**.

For optimal performance, the pipeline is designed to run on a **Linux HPC system** with a **SLURM scheduler**, although it can also run on a multi-core local machine with Singularity.  

Disk space requirements depend on input data size; for typical analyses with Atlantic salmon datasets and 12 samples, at least **200 GB of free storage** is recommended.
  
Memory and CPU usage may scale with the number of samples and genomic data size.


## 3. Input Arguments

| Argument        | Description                     | Required                    |
| --------------- | ------------------------------- | --------------------------- |
| `--genome`      | Reference genome FASTA          | Yes                         |
| `--annotation`  | Gene annotation GTF             | Yes                         |
| `--atac_bam`    | ATAC-seq BAM files              | Optional (for ATAC-QTL)    |
| `--rna_bam`     | RNA-seq BAM files               | Optional (for eQTL)        |
| `--genotype`    | Genotype VCF                    | Yes                         |
| `--meta`        | Sample metadata CSV             | Yes                         |
| `--outdir`      | Output directory                | Yes                         |
| `--ld_genotype` | Optional VCF for LD analysis    | No                          |
| `--fpkm_cutoff` | Filter threshold for expression | Optional                    |
| `--atac_qtl`    | Enable ATAC-QTL analysis        | Yes (true/false)            |
| `--eqtl_qtl`    | Enable eQTL analysis            | Yes (true/false)            |
| `--external_ld` | Use external LD genotype file   | Yes (true/false)            |

---

## 4. Data Input Explanation

| Data Type            | Description                     | Example path                             |
| -------------------- | ------------------------------- | ---------------------------------------- |
| **Reference genome** | FASTA file of the genome        | `data/ref/genome.fa`                     |
| **Annotation**       | Gene annotation GTF             | `data/ref/annotation.gtf`                |
| **ATAC BAM**         | Aligned ATAC-seq reads          | `data/atac_bam/*.bam`                    |
| **RNA BAM**          | Aligned RNA-seq reads           | `data/rna_bam/*.bam`                     |
| **ATAC count**       | Consensus peak counts           | `data/atac_consensus_peak_featureCounts.txt` |
| **RNA count**        | Gene-level counts               | `data/rna_gene_level_count_salmon.txt`   |
| **Genotype VCF**     | Genotype calls for QTL mapping  | `data/genotype.vcf.gz`                   |
| **Metadata CSV**     | Sample information              | `data/meta/Brain.csv`                     |
| **LD genotype VCF**  | Optional, for LD-based analyses | `data/ld_genotype.vcf.gz`                |

---   



## 5. Example Running

### Run with ATAC and RNA-seq on HPC

Note: you may need to customize the nextflow.config file to make it works with your HPC.

```bash
nextflow run nf-rasqual/main.nf -profile hpc \
  --genome data/ref/genome.fa \
  --annotation data/ref/annotation.gtf \
  --atac_bam data/atac_bam/*.bam \
  --rna_bam data/rna_bam/*.bam \
  --genotype data/genotype.vcf.gz \
  --meta data/meta/Brain.csv \
  --outdir results \
  --atac_qtl true \
  --eqtl_qtl true \
  --external_ld false ## not using external LD

```

### Run with only RNA-seq on HPC

Note: you may need to customize the nextflow.config file to make it works with your HPC.

```bash
nextflow run nf-rasqual/main.nf -profile hpc \
  --genome data/ref/genome.fa \
  --annotation data/ref/annotation.gtf \
  --rna_bam data/rna_bam/*.bam \
  --genotype data/genotype.vcf.gz \
  --meta data/meta/Brain.csv \
  --outdir results \
  --atac_qtl false \
  --eqtl_qtl true \
  --external_ld false ## not using external LD
```


### Notes

* The workflow automatically splits VCFs and counts per chromosome for parallel RASQUAL runs.
* Users can specify FPKM cutoffs and other parameters as needed.
* Multi-core CPUs accelerate RASQUAL runs; GPUs are not required.

---

If you want, I can also add a **supplementary table mapping each workflow step to memory, CPU, tool, input/output, and notes** directly under the README so it fully addresses reviewer reproducibility concerns. This would make it very easy to follow and reproduce.

Do you want me to do that next?

## Running the Pipeline for Atlantic Salmon

## Download Reference and Annotation Files

All analyses use the Atlantic salmon genome assembly **Ssal v3.1** and **Ensembl annotation release 106**. You can download the required files from the Ensembl FTP site:

- Genome assembly (FASTA): [https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz](https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/)
- Annotation (GTF): [https://ftp.ensembl.org/pub/release-106/gtf/salmo_salar/Salmo_salar.Ssal_v3.1.106.gtf.gz](https://ftp.ensembl.org/pub/release-106/gtf/salmo_salar/)

Ensure that both files are saved locally and paths are correctly set for downstream analyses.

To run the pipeline, configure the input files and parameters in the pipeline script, then execute the following scripts:

- run_officical_brain.sh
- run_officical_liver.sh
- run_officical_gonad.sh
- run_officical_muscle.sh
- run_officical_gill.sh
