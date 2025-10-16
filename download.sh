#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<EOF
This script downloads the bam files and feature count files used in the nf-rasqual paper from salmobase.org.

Usage: $(basename "$0") -t <Tissue> -a <Assay> [--counts-only] [--check]

Required:
  -t, --tissue   Tissue to download (Brain, Gill, Gonad, Liver, Muscle)
  -a, --assay    Assay to download (ATAC or RNA)

Optional:
  --counts-only  Download only the count matrix (skip BAM/BAM index files)
  --check        Dry run; only verify that URLs are reachable
  -h, --help     Show this help message
EOF
  exit "${1:-0}"
}

valid_tissues=(Brain Gill Gonad Liver Muscle)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="${script_dir}/data"
meta_base="${data_dir}/meta"

tissue=""
assay=""
check_only=false
counts_only=false

tolower() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

toupper() {
  printf '%s' "$1" | tr '[:lower:]' '[:upper:]'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--tissue)
      [[ $# -lt 2 ]] && { echo "Error: Missing value for $1." >&2; usage 1; }
      tissue="$2"
      shift 2
      ;;
    -a|--assay)
      [[ $# -lt 2 ]] && { echo "Error: Missing value for $1." >&2; usage 1; }
      assay="$2"
      shift 2
      ;;
    --check)
      check_only=true
      shift
      ;;
    -c|--counts-only)
      counts_only=true
      shift
      ;;
    -h|--help)
      usage 0
      ;;
    *)
      echo "Error: Unknown argument: $1" >&2
      usage 1
      ;;
  esac
done

if [[ -z "${tissue}" || -z "${assay}" ]]; then
  echo "Error: Both tissue and assay must be specified." >&2
  usage 1
fi

# Normalize tissue input and validate
matched_tissue=""
tissue_lower="$(tolower "${tissue}")"
for option in "${valid_tissues[@]}"; do
  if [[ "${tissue_lower}" == "$(tolower "${option}")" ]]; then
    matched_tissue="${option}"
    break
  fi
done

if [[ -z "${matched_tissue}" ]]; then
  echo "Error: Unsupported tissue '${tissue}'. Valid options: ${valid_tissues[*]}." >&2
  exit 1
fi
tissue="${matched_tissue}"

assay_upper="$(toupper "${assay}")"
if [[ "${assay_upper}" != "ATAC" && "${assay_upper}" != "RNA" ]]; then
  echo "Error: Unsupported assay '${assay}'. Valid options: ATAC, RNA." >&2
  exit 1
fi

meta_file="${meta_base}/${tissue}.csv"
if [[ "${counts_only}" == false && ! -f "${meta_file}" ]]; then
  echo "Error: Meta file not found for tissue '${tissue}' at ${meta_file}." >&2
  exit 1
fi

fetch() {
  local url=$1
  local destination=$2

  if "${check_only}"; then
    echo "Checking ${url}"
    if ! curl -kfsI "${url}" -o /dev/null; then
      echo "Link check failed for ${url}" >&2
      exit 1
    fi
    return
  fi

  mkdir -p "$(dirname "${destination}")"
  echo "Downloading ${url} -> ${destination}"
  if ! curl -k -fL "${url}" -o "${destination}"; then
    echo "Failed to download ${url}" >&2
    exit 1
  fi
}

case "${assay_upper}" in
  ATAC)
    atac_base="https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/ATAC"
    count_url="${atac_base}/${tissue}/results/bwa/mergedLibrary/macs/narrowPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt"
    fetch "${count_url}" "${data_dir}/${tissue}/atac_consensus_peak_counts.txt"

    if [[ "${counts_only}" == false ]]; then
      tail -n +2 "${meta_file}" | while IFS=, read -r genotype_id sample_description atac_count_id atac_bam_id rna_count_id rna_bam_id sex stage; do
        atac_bam_id="${atac_bam_id//$'\r'/}"
        [[ -z "${atac_bam_id}" ]] && continue
        bam_url="${atac_base}/${tissue}/results/bwa/mergedLibrary/${atac_bam_id}"
        fetch "${bam_url}" "${data_dir}/${tissue}/atac_bam/${atac_bam_id}"
        bai_url="${bam_url}.bai"
        fetch "${bai_url}" "${data_dir}/${tissue}/atac_bam/${atac_bam_id}.bai"
      done
    fi
    ;;
  RNA)
    rna_base="https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/RNA"
    count_url="${rna_base}/${tissue}/results/star_salmon/salmon.merged.gene_counts.tsv"
    fetch "${count_url}" "${data_dir}/${tissue}/rna_gene_counts.tsv"

    if [[ "${counts_only}" == false ]]; then
      tail -n +2 "${meta_file}" | while IFS=, read -r genotype_id sample_description atac_count_id atac_bam_id rna_count_id rna_bam_id sex stage; do
        rna_bam_id="${rna_bam_id//$'\r'/}"
        [[ -z "${rna_bam_id}" ]] && continue
        bam_url="${rna_base}/${tissue}/results/star_rsem/${rna_bam_id}"
        fetch "${bam_url}" "${data_dir}/${tissue}/rna_bam/${rna_bam_id}"
        bai_url="${bam_url}.bai"
        fetch "${bai_url}" "${data_dir}/${tissue}/rna_bam/${rna_bam_id}.bai"
      done
    fi
    ;;
esac
