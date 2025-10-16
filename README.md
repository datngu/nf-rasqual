# nf-rasqual Pipeline

## Overview

This pipeline processes and analyzes genomic, transcriptomic, and ATAC-seq data for QTL (Quantitative Trait Loci) analysis. It uses genotype and phenotype data to identify associations, focusing on eQTLs (expression QTLs) and ATAC-QTLs. The pipeline is designed for flexibility, allowing users to configure various parameters, including the use of external linkage disequilibrium (LD) data for multiple testing correction with EigenMT.



## 2. Installation and Dependencies

### Dependencies

* **Nextflow** (v24.4 or higher)  
* **Singularity** for containerized execution  

Since pipeline is implemented with Nextflow (DSL2), you will need Nextflow to run it. Installation instructions are available on the Nextflow homepage: [https://www.nextflow.io/](https://www.nextflow.io/).  

All tool dependencies are included in the corresponding container (next section), so running with Docker or Singularity is sufficient.  

By default, nf-rasqual can run with Singularity on an HPC using SLURM (**-profile saga,singularity**). You may need to adjust **nextflow.config** to match your system configuration.


Please keep in mind that the provided `nextflow.config` is a template and require modifications to suit your specific HPC environment. Adjust parameters such as **executor settings, queue names, account, etc** according to your system's requirements.

For further customization, consult the Nextflow documentation or your system admin, we can't prodive instruction on this: [https://www.nextflow.io/docs/latest/config.html](https://www.nextflow.io/docs/latest/config.html).


### Container

A pre-built Docker container includes all required software (RASQUAL, R scripts, Python scripts, bcftools, etc.) that will be automaticly downloaded  by **Nextflow**

```
ndatth/rasqual:v0.0.0
```

Alternatively, if interest, users can build the container themselves as following the instructions and using the source code and libraries provided in the `container` directory of the repository: [https://github.com/datngu/nf-rasqual/tree/main/container](https://github.com/datngu/nf-rasqual/tree/main/container).


### Hardware Requirements

The most resource-intensive step of the pipeline requires **64 GB of RAM** and **16 CPU cores**.

For optimal performance, the pipeline is designed to run on a **Linux HPC system** with a **SLURM scheduler**, although it can also run on a multi-core local machine with Singularity (not recommended and not tested yet).

Disk space requirements depend on input data size; for typical analyses with Atlantic salmon datasets and 12 samples, at least **500 GB of free storage** is recommended.
  
Memory and CPU usage may scale with the number of samples and genomic data size.


## Input files and parameters

| Parameter         | Required                         | Description                                                       | Example/default value                          |
| ----------------- | -------------------------------- | ----------------------------------------------------------------- | -------------------------------------- |
| `--genome`        | Yes                              | Reference genome FASTA used when preparing rasqual inputs.        | `data/ref/genome.fa`                   |
| `--annotation`    | Yes                              | Gene annotation GTF parsed for RNA gene definitions.              | `data/ref/annotation.gtf`              |
| `--atac_bam`      | Conditional (`--atac_qtl true`)  | ATAC BAM files used for allele-specific counting.                 | `data/atac_bam/*.bam`                  |
| `--atac_count`    | Conditional (`--atac_qtl true`)  | ATAC peak count matrix filtered before rasqual.                   | `data/atac_consensus_peak_featureCounts.txt` |
| `--rna_bam`       | Conditional (`--eqtl_qtl true`)  | RNA BAM files used for allele-specific counting.                  | `data/rna_bam/*.bam`                   |
| `--rna_count`     | Conditional (`--eqtl_qtl true`)  | RNA gene-level counts filtered before rasqual.                    | `data/rna_gene_level_count_salmon.txt` |
| `--genotype`      | Yes                              | Sample genotype VCF driving cis-QTL discovery.                    | `data/genotype.vcf.gz`                 |
| `--ld_genotype`   | Conditional (`--external_ld true`)  | External LD panel if running eigenMT with external LD.            | `data/ld_genotype.vcf.gz`              |
| `--meta`          | Yes                              | Sample metadata linking assays and covariates.                    | `data/meta/Brain.csv`                  |
| `--outdir`        | Optional                         | Final output directory for merged results.                        | `results` (default)                             |
| `--trace_dir`     | Optional                         | Directory capturing intermediate `publishDir` artefacts.          | `trace_dir` (default)                           |
| `--chrom`         | Optional                         | Chromosome identifiers to process.                                | `1..29`(default)                                |
| `--phenotype_PCs` | Optional                         | Number of phenotype PCs to include as covariates.                 | `2`(default)                                    |
| `--exp_prop`      | Optional                         | Fraction of samples that must pass the FPKM cutoff.               | `0.5` (default)                                  |
| `--fpkm_cutoff`   | Optional                         | Minimum FPKM threshold applied during expression filtering.       | `0.5` (default)                                 |
| `--atac_window`   | Optional               | Cis distance (bp) for ATAC rasqual and eigenMT pairing.           | `10000` (default)                               |
| `--eqtl_window`   | Optional                   | Cis distance (bp) for RNA rasqual and eigenMT pairing.            | `500000` (default)                              |
| `--atac_qtl`      | Optional                         | Toggle to run the ATAC-QTL branch. `true` or `false`                | `true` (default)                                  |
| `--eqtl_qtl`      | Optional                         | Toggle to run the eQTL branch. `true` or `false`                                   | `true` (default)                                  |
| `--external_ld`   | Optional                         | Enable external LD-aware eigenMT processing. `true` or `false`                     | `false` (default)                                |

**BAM indexes:** Every BAM provided to the workflow must have a matching `.bam.bai` index in the same directory; the pipeline copies these alongside the BAMs.

**Chromosome identifiers:** Supply `--chrom` values exactly as they appear in your count tables and VCF (for example `chr1`). For named contigs, pass a Groovy list on the command line (`--chrom '["chr1","chr2"]'`) or set `params.chrom = ["chr1","chr2"]` in a custom `nextflow.config`.

### Sample Metadata CSV Format

The metadata CSV links each sample across genotype, RNA-seq, and ATAC-seq inputs and provides covariates used during QTL modeling. The file must include a header row with the following columns:

| Column               | Required    | Description                                                                                   | Example value                                                         |
| -------------------- | ----------- | --------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `genotype_id`        | Yes         | Sample identifier exactly matching the column name in the genotype VCF.                       | `A14_AF`                                                              |
| `sample_description` | Yes         | Human-readable label used in logs and reports.                                                | `MatureFemale_R1`                                                     |
| `atac_count_id`      | Conditional | Column name for the sample in the ATAC count matrix when `--atac_qtl true`; leave empty otherwise.| `AtanticSalmon_ATAC_MatureFemale_Brain_R1.mLb.clN.bam`                |
| `atac_bam_id`        | Conditional | Filename of the ATAC BAM when `--atac_qtl true`; leave empty otherwise.                      | `AtlanticSalmon_ATAC_Brain_Mature_Female_R1.mLb.clN.sorted.bam`       |
| `rna_count_id`       | Conditional | Column name for the sample in the RNA count matrix when `--eqtl_qtl true`; leave empty otherwise.    | `AtlanticSalmon_RNA_Brain_Mature_Female_R1`                           |
| `rna_bam_id`         | Conditional | Filename of the RNA BAM when `--eqtl_qtl true`; leave empty otherwise.                        | `AtlanticSalmon_RNA_Brain_Mature_Female_R1.markdup.sorted.bam`        |
| `sex`                | Optional    | Numeric covariate passed to RASQUAL (e.g. `0` = female, `1` = male).                           | `0`                                                                   |
| `stage`              | Optional    | Additional covariate (e.g. `0` = immature, `1` = mature); add extra covariate columns as needed.| `1`                                                                   |

Each row therefore ties the genotype column (e.g. `A14_AF`) to the matching BAM and count identifiers representing the same individual. When rasqual covariates are constructed the pipeline keeps those first six identifier columns fixed and appends every additional column from the metadata (for example `sex`, `stage`, or any custom numeric covariate you add). To introduce extra covariates simply add extra columns to the CSV.

### ATAC feature count format

nf-core/atacseq produces a consensus peak feature count table (featureCounts output) that matches the expected ATAC input for this workflow. Files are tab-separated, begin with optional comment lines starting with `#`, and include a header with six feature columns (`Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`) followed by one column per ATAC sample. Peak rows contain integer counts for each sample. See example `data/atac_consensus_peak_featureCounts.txt`.

### RNA feature count format

nf-core/rnaseq emits merged gene-level raw read count matrices that serve as the RNA input here. These tables are tab-separated with a header containing `gene_id`, `gene_name`, and one column per RNA sample. Values are typically integer counts of the number of reads mapping to the gene. See example `data/rna_gene_level_count_salmon.txt`.

## 5. A copy-paste example to running the pipeline

To demostrate the pipeline ultility, we provide a detailed copy-paste tutorial to test/reproduce the running for Brain tissues of for Atlantic Salmon.

### Download Reference and Annotation Files

All analyses use the Atlantic salmon genome assembly **Ssal v3.1** and **Ensembl annotation release 106**. You can download the required files from the Ensembl FTP site:

- Genome assembly (FASTA): [https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz](https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/)
- Annotation (GTF): [https://ftp.ensembl.org/pub/release-106/gtf/salmo_salar/Salmo_salar.Ssal_v3.1.106.gtf.gz](https://ftp.ensembl.org/pub/release-106/gtf/salmo_salar/)

Ensure that both files are saved locally and paths are correctly set for downstream analyses. The snippet below downloads and prepares the reference in `data/ref/`:

```bash

# git clone git@github.com:datngu/nf-rasqual.git
# cd nf-rasqual


mkdir -p data/ref
curl -Lk https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz \
  -o data/ref/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz
curl -Lk https://ftp.ensembl.org/pub/release-106/gtf/salmo_salar/Salmo_salar.Ssal_v3.1.106.gtf.gz \
  -o data/ref/Salmo_salar.Ssal_v3.1.106.gtf.gz
gunzip -k data/ref/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz
gunzip -k data/ref/Salmo_salar.Ssal_v3.1.106.gtf.gz

```

### Download LD Genotype Panel

The nf-rasqual paper uses an external LD genotype panel hosted on SalmoBase. Download it to `data/ref/`:

```bash
mkdir -p data/ref
curl -Lk https://salmobase.org/datafiles/datasets/nf-rasqual-data/ALL_chrome_phased_filtered_HWE_1e6_biSNPs_MAF_0.01.vcf.gz \
  -o data/ref/ALL_chrome_phased_filtered_HWE_1e6_biSNPs_MAF_0.01.vcf.gz

```

## Download BAM and Count Files

The repository ships with a helper script that fetches the BAM (and .bai) and feature-count files from Salmobase.org for the nf-rasqual Atlantic salmon analysis. From the repository root you can download the Brain dataset (ATAC or RNA) as follows:

```bash
## the download may take a few hours depending on your internet connection
# Download ATAC files (consensus peak counts + BAMs)
bash download.sh -t Brain -a ATAC

## the download may take a few hours depending on your internet connection
# Download RNA-seq files (Salmon counts + BAMs)
bash download.sh -t Brain -a RNA

```

Downloaded artefacts are written under `data/Brain/` in assay-specific subdirectories.

### Run the Brain Analysis on HPC

Once the reference, BAMs, and count files are in place, you can launch the Brain ATAC/RNA run on a SLURM-based HPC using Nextflow. Adjust module loads, work directories, and resource requirements as needed for your environment.

```bash
## on Saga HPC https://documentation.sigma2.no/hpc_machines/saga.html#saga

module load Nextflow/24.04.2
## Singularity is usually available by default on Saga
# If needed, load a specific Singularity version:
# module load singularity/your_version or contact your sysadmin :D

export NXF_SINGULARITY_CACHEDIR=$PWD/container

nextflow run main.nf -profile saga,singularity \
  --meta "$PWD/data/meta/Brain.csv" \
  --genome "$PWD/data/ref/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa" \
  --annotation "$PWD/data/ref/Salmo_salar.Ssal_v3.1.106.gtf" \
  --atac_bam "$PWD/data/Brain/atac_bam/*.bam" \
  --atac_count "$PWD/data/Brain/atac_consensus_peak_counts.txt" \
  --rna_bam "$PWD/data/Brain/rna_bam/*.bam" \
  --rna_count "$PWD/data/Brain/rna_gene_counts.tsv" \
  --genotype "$PWD/data/genotype.vcf.gz" \
  --outdir "$PWD/results/Brain" \
  --atac_window 10000 \
  --external_ld true \
  --ld_genotype "$PWD/data/ref/ALL_chrome_phased_filtered_HWE_1e6_biSNPs_MAF_0.01.vcf.gz" \
  --atac_qtl true \
  --eqtl_qtl true

```

**FOR REVIEWERS:** If you are testing/running the pipeline on your system. The provided `nextflow.config` is a template and may require modifications to suit your specific HPC environment. Adjust parameters such as **executor settings, queue names, account, etc** according to your system's requirements. We can't help with your specific setup, but we recommend consulting the Nextflow documentation or your system admin for further customization: [https://www.nextflow.io/docs/latest/config.html](https://www.nextflow.io/docs/latest/config.html).


The example above runs both ATAC-QTL and eQTL components. Remove `--atac_bam`, `--atac_count`, and `--atac_qtl true` if you only wish to execute the RNA component. Repeat the same download and execution steps for the other tissues (Gill, Gonad, Liver, Muscle) by swapping the tissue name in the commands; doing so reproduces the analyses performed in the nf-rasqual paper. Adjust `--outdir`, and `-w` (Nextflow work directory) as needed to point to project-specific locations with sufficient storage.
