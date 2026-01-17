#!/bin/bash

set -euo pipefail

# === Usage check ===
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <species_name> <reference_name> <base_output_dir> [--no-cache] [--mapq <value>] [--remove-temp] [--aligner bwa|minimap2|bbmap] [--aligner-cmd <command>]"
    exit 1
fi

# === Required arguments ===
SPECIES=$1
REFERENCE=$2
BASE_DIR=$3
shift 3

# === Optional flags ===
NO_CACHE=false
REMOVE_TEMP=false
MAPQ=""
ALIGNER=""
ALIGNER_CMD=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-cache)
            NO_CACHE=true
            shift
            ;;
        --remove-temp)
            REMOVE_TEMP=true
            shift
            ;;
        --mapq)
            MAPQ="$2"
            shift 2
            ;;
        --aligner)
            ALIGNER="$2"
            shift 2
            ;;
        --aligner-cmd)
            ALIGNER_CMD="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# === Set default to bwa if neither provided ===
if [[ -z "$ALIGNER_CMD" && -z "$ALIGNER" ]]; then
    echo "No aligner specified. Defaulting to 'bwa'"
    ALIGNER="bwa"
fi

# === Derive ALIGNER_CMD from ALIGNER if needed ===
if [[ -z "$ALIGNER_CMD" && -n "$ALIGNER" ]]; then
    case "$ALIGNER" in
        bwa)
            ALIGNER_CMD="bwa mem -t {cores} {ref} {fq}"
            ;;
        minimap2)
            ALIGNER_CMD="minimap2 -t {cores} -ax sr {ref} {fq}"
            ;;
        bbmap)
            ALIGNER_CMD="bbmap.sh ref={ref} threads={cores} in={fq} out=stdout.sam"
            ;;
        *)
            echo "Unsupported aligner: $ALIGNER"
            exit 1
            ;;
    esac
fi

# === Final validation for placeholders ===
if [[ "$ALIGNER_CMD" != *"{ref}"* || "$ALIGNER_CMD" != *"{fq}"* || "$ALIGNER_CMD" != *"{cores}"* ]]; then
    echo "$ALIGNER_CMD"
    echo "--aligner-cmd must include placeholders: {ref}, {fq}, and {cores}"
    exit 1
fi

# === Paths ===
SPECIES_FASTA="${BASE_DIR}/${SPECIES}/${SPECIES}.fasta"
REFERENCE_FASTA="${BASE_DIR}/${REFERENCE}/${REFERENCE}.fasta"
BAM_DIR="${BASE_DIR}/BAMs"
PLOTS_DIR="${BASE_DIR}/Plots"
FASTQ="${BASE_DIR}/${SPECIES}/${SPECIES}_reads_to_${REFERENCE}.fastq"
RAW_BAM="${BAM_DIR}/${SPECIES}_to_${REFERENCE}_raw.bam"
BAM="${BAM_DIR}/${SPECIES}_to_${REFERENCE}.bam"
LOG="${BAM%.bam}.log"

# === Input checks ===
if [[ ! -f "$SPECIES_FASTA" ]]; then
    echo "Species FASTA not found: $SPECIES_FASTA"
    exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "Reference FASTA not found: $REFERENCE_FASTA"
    exit 1
fi

# === Final BAM caching ===
if [[ -f "$BAM" && "$NO_CACHE" == false ]]; then
    echo "Final BAM already exists: $BAM"
    exit 0
fi

mkdir -p "$BAM_DIR"

# === FASTQ generation ===
if [[ -f "$FASTQ" && "$NO_CACHE" == false ]]; then
    echo "Reusing existing FASTQ file: $FASTQ"
else
    echo "Generating FASTQ from $SPECIES_FASTA to $FASTQ"
    python create_fastq_fragments.py "$SPECIES_FASTA" --output "$FASTQ" --force
fi

# === Prepare filter command ===
FILTER_CMD="python filter_sam.py -"
if [[ -n "$MAPQ" ]]; then
    FILTER_CMD+=" --mapq $MAPQ"
    FILTER_CMD+=" --mapq-hist-folder $PLOTS_DIR"
fi

CORES=$(nproc)

# === Substitute aligner placeholders ===
ALIGN_CMD=${ALIGNER_CMD//\{ref\}/$REFERENCE_FASTA}
ALIGN_CMD=${ALIGN_CMD//\{fq\}/$FASTQ}
ALIGN_CMD=${ALIGN_CMD//\{cores\}/$CORES}

# === Alignment and filtering logic ===
if [[ "$REMOVE_TEMP" == false ]]; then
    if [[ -f "$RAW_BAM" && "$NO_CACHE" == false ]]; then
        echo "Using cached raw BAM: $RAW_BAM"
    else
        echo "Aligning and saving raw BAM..."
        eval "$ALIGN_CMD" | samtools view -h -o "$RAW_BAM"
        echo "Raw BAM saved: $RAW_BAM"
    fi

    echo "Filtering $RAW_BAM to produce final BAM..."
    samtools view -h "$RAW_BAM" \
        | eval "$FILTER_CMD" 2> "$LOG" \
        | samtools sort -@ "$CORES" -o "$BAM"
else
    echo "Aligning, filtering, and sorting in stream..."
    eval "$ALIGN_CMD" \
        | eval "$FILTER_CMD" 2> "$LOG" \
        | samtools sort -@ "$CORES" -o "$BAM"

    [[ -f "$FASTQ" ]] && rm "$FASTQ"
fi

# === Index final BAM ===
echo "Indexing final BAM..."
samtools index "$BAM"

# === Done ===
echo "Finished: $BAM"
echo "Filter stats written to: $LOG"
[[ "$REMOVE_TEMP" == false ]] && echo "Unfiltered alignment saved: $RAW_BAM"
