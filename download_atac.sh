#!/usr/bin/env bash
set -euo pipefail

tis="$1"
CURRENT_DIR="$(pwd)"

# Create output directories
mkdir -p "${CURRENT_DIR}/data/${tis}/atac_bam"

# Download count matrix
wget --no-check-certificate -q \
  "https://salmobase.org/datafiles/datasets/nf-rasqual-data/${tis}/atac_consensus_peak_featureCounts.txt" \
  -O "${CURRENT_DIR}/data/${tis}/atac_consensus_peak_featureCounts.txt"

# Download ATAC BAM and BAI files
wget -r -np -nH -R "index.html*" --cut-dirs=5 -e robots=off --no-check-certificate -q \
  -P "${CURRENT_DIR}/data/${tis}/atac_bam" \
  "https://salmobase.org/datafiles/datasets/nf-rasqual-data/${tis}/atac_bam/"
