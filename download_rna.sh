#!/usr/bin/env bash
set -euo pipefail

tis="$1"
CURRENT_DIR="$(pwd)"

# Create output directories
mkdir -p "${CURRENT_DIR}/data/${tis}/rna_bam"

# Download count matrix
wget --no-check-certificate -q \
  "https://salmobase.org/datafiles/datasets/nf-rasqual-data/${tis}/rna_gene_level_count_salmon.txt" \
  -O "${CURRENT_DIR}/data/${tis}/rna_gene_level_count_salmon.txt"

# Download RNA BAM and BAI files
wget -r -np -nH -R "index.html*" --cut-dirs=5 -e robots=off --no-check-certificate -q \
  -P "${CURRENT_DIR}/data/${tis}/rna_bam" \
  "https://salmobase.org/datafiles/datasets/nf-rasqual-data/${tis}/rna_bam/"
