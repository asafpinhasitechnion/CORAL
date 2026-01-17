#!/bin/bash

set -euo pipefail

# === Usage ===
if [[ "$#" -lt 2 ]]; then
    echo "Usage: $0 <species_name> <base_output_dir> [--no-cache]"
    exit 1
fi

# === Input Arguments ===
SPECIES_NAME=$1
BASE_DIR=$2
NO_CACHE=${3:-""}
FASTA_FILE="${BASE_DIR}/${SPECIES_NAME}/${SPECIES_NAME}.fasta"

# === Check Tools ===
command -v bwa >/dev/null 2>&1 || { echo >&2 "bwa not found in PATH."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools not found in PATH."; exit 1; }

# === Check FASTA ===
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "FASTA file not found: $FASTA_FILE"
    exit 1
fi

# === Check for Existing Index Files ===
BWA_INDEX_DONE=true
SAMTOOLS_INDEX_DONE=true

for ext in amb ann bwt pac sa; do
    [[ ! -f "${FASTA_FILE}.${ext}" ]] && BWA_INDEX_DONE=false
done

[[ ! -f "${FASTA_FILE}.fai" ]] && SAMTOOLS_INDEX_DONE=false

if [[ "$NO_CACHE" != "--no-cache" && "$BWA_INDEX_DONE" == true && "$SAMTOOLS_INDEX_DONE" == true ]]; then
    echo "Index files already exist for $SPECIES_NAME. Skipping indexing."
    exit 0
fi

# === Clean Old Indices ===
rm -f "${FASTA_FILE}".{amb,ann,bwt,pac,sa,fai}

# === Run Indexing ===
echo "Indexing with BWA..."
bwa index "$FASTA_FILE" || { echo "BWA indexing failed."; exit 1; }

echo "Indexing with Samtools..."
samtools faidx "$FASTA_FILE" || { echo "Samtools indexing failed."; exit 1; }

echo "Indexing complete for $SPECIES_NAME in $BASE_DIR"
