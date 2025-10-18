#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: ${0##*/} TISSUE

Download ATAC count matrix and BAM files for a given tissue from salmobase.

Positional arguments:
  TISSUE    Tissue name (e.g., Brain, Gill, Gonad, Liver, Muscle)

Example:
  ${0##*/} Brain

The downloaded data will be saved in the following directory structure:
  ./data/TISSUE/atac_bam/*.{bam,.bai}
  ./data/TISSUE/atac_consensus_peak_featureCounts.txt
EOF
}

if [[ $# -ne 1 ]]; then
  usage >&2
  exit 1
fi

case "$1" in
  -h|--help)
    usage
    exit 0
    ;;
esac

tis="$1"
CURRENT_DIR="$(pwd)"

# Create output directories
mkdir -p "${CURRENT_DIR}/data/${tis}/atac_bam"

# Download count matrix
wget --no-check-certificate \
  "https://salmobase.org/datafiles/datasets/nf-rasqual-data/${tis}/atac_consensus_peak_featureCounts.txt" \
  -O "${CURRENT_DIR}/data/${tis}/atac_consensus_peak_featureCounts.txt"

# Download ATAC BAM and BAI files
wget -r -np -nH -R "index.html*" --cut-dirs=5 -e robots=off --no-check-certificate \
  -P "${CURRENT_DIR}/data/${tis}/atac_bam" \
  "https://salmobase.org/datafiles/datasets/nf-rasqual-data/${tis}/atac_bam/"
