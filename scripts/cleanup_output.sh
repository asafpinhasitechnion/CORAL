#!/bin/bash

set -euo pipefail

# === USAGE ===
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <base_output_dir> [--bams] [--pileup] [--intervals] [--genomes] --t1 <taxa1> --t2 <taxa2> --out <outgroup>"
  exit 1
fi

BASE_DIR=$1
shift

# === FLAGS ===
DELETE_BAMS=false
DELETE_PILEUP=false
DELETE_INTERVALS=false
DELETE_GENOMES=false

T1_NAME=""
T2_NAME=""
OUT_NAME=""

# === PARSE ARGUMENTS ===
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bams) DELETE_BAMS=true ; shift ;;
    --pileup) DELETE_PILEUP=true ; shift ;;
    --intervals) DELETE_INTERVALS=true ; shift ;;
    --genomes) DELETE_GENOMES=true ; shift ;;
    --t1) T1_NAME="$2" ; shift 2 ;;
    --t2) T2_NAME="$2" ; shift 2 ;;
    --out) OUT_NAME="$2" ; shift 2 ;;
    *) echo "Unknown option: $1" ; exit 1 ;;
  esac
done

# === DELETE BAM FILES ===
if $DELETE_BAMS; then
  echo "Removing BAM files..."
  rm -rf "$BASE_DIR/BAMs"
fi

# === DELETE PILEUP FILE ===
if $DELETE_PILEUP; then
  echo "Removing pileup file..."
  rm -f "$BASE_DIR/${OUT_NAME}__${T1_NAME}__${T2_NAME}.pileup.gz"
fi

# === DELETE INTERVALS ===
if $DELETE_INTERVALS; then
  echo "Removing interval coverage files..."
  rm -rf "$BASE_DIR/Intervals"
fi

# === DELETE GENOMES ===
if $DELETE_GENOMES; then
  echo "Removing genome FASTA and index files..."
  for genome in "$BASE_DIR/$T1_NAME" "$BASE_DIR/$T2_NAME" "$BASE_DIR/$OUT_NAME"; do
    if [[ -d "$genome" ]]; then
      echo "   - Cleaning $genome"
      rm -rf "$genome"
    fi
  done
fi

echo "Cleanup complete."
