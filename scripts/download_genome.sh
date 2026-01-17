#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <genome_name> <genome_accession> <output_directory> [--no-cache]"
    exit 1
fi

GENOME_NAME=$1
ACCESSION=$2
BASE_OUTPUT_DIR=$3
NO_CACHE=${4:-false}

OUTPUT_DIR="${BASE_OUTPUT_DIR}/${GENOME_NAME}"
FASTA_FILE="${OUTPUT_DIR}/${GENOME_NAME}.fasta"

# === Check for existing genome ===
if [[ -f "$FASTA_FILE" && "$NO_CACHE" != "--no-cache" ]]; then
    echo "Cached genome found for $GENOME_NAME. Skipping download."
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
TEMP_DIR="${BASE_OUTPUT_DIR}/temp_${GENOME_NAME}"
mkdir -p "$TEMP_DIR"

echo "Downloading genome data for accession: $ACCESSION"
datasets download genome accession "$ACCESSION" --filename "${TEMP_DIR}/${GENOME_NAME}.zip"

echo "Unzipping genome..."
unzip -q "${TEMP_DIR}/${GENOME_NAME}.zip" -d "$TEMP_DIR"

FASTA_SOURCE=$(find "$TEMP_DIR" -name "*.fna" -print -quit)

if [ -z "$FASTA_SOURCE" ]; then
    echo "No FASTA file found in the downloaded data."
    rm -rf "$TEMP_DIR"
    exit 1
fi

mv "$FASTA_SOURCE" "$FASTA_FILE"
echo "FASTA saved at: $FASTA_FILE"

echo "Stripping FASTA headers (keeping only first word)..."
sed -i 's/ .*$//' "$FASTA_FILE"

rm -rf "$TEMP_DIR"
