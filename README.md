# CORAL:  Comparative Orthologous Read-based Analysis of Lineage Substitutions

CORAL is a tool for scalable extraction, detection, and analysis of point mutations across species evolutionary history.
It aligns multiple species to a shared reference genome, simulates reads, filters alignments by mapping quality, extracts unambiguous trinucleotide substitutions, and summarizes mutation rates and spectra.

Reference  
Please cite: will be available upon publication

---
<img width="8131" height="2352" alt="Figure1_Species_mutations_Asaf (8)" src="https://github.com/user-attachments/assets/80d9c5af-6981-4d76-8795-cb0a9b882027" />

## Core Concepts

### Mutation Context Extraction
CORAL identifies genomic positions where:
- All species share the same flanking bases
- Exactly one lineage differs at the central base

These sites represent unambiguous point mutations in trinucleotide context.
They form the basis for mutation spectrum analysis, mutation rate estimation, and genome-wide summaries.

---

### Distributed Alignment Using Pseudo-Reads
Each non-outgroup species is aligned independently to a shared reference (outgroup) genome using simulated short reads.

This approach:
- Avoids full multiple sequence alignment
- Minimizes reference bias associated with MSAs
- Is computationally lightweight
- Scales naturally across many species

Standard short-read aligners (for example, BWA) are used for efficiency and robustness.

---

## Functional Workflow

### Step 1: Genome Preparation
- Download genomes by NCBI assembly accession
- Index the reference genome for alignment

### Step 2: Read Simulation and Alignment
- Simulate FASTQ reads by sliding a window across genomes
- Align simulated reads to the outgroup reference
- Filter alignments by MAPQ and coverage
- Allow customization of aligner and parameters

### Step 3: Mutation Detection
- Generate pileups from reference and aligned BAMs
- Extract unambiguous trinucleotide substitutions
- Optionally retain genomic positions

### Step 4: Normalization and Analysis
- Normalize mutation counts by underlying trinucleotide abundance
- Collapse complementary strands into canonical spectra
- Generate summary tables and visualizations

---

## Example Usage

### Single-lineage analysis (outgroup + two ingroups)
coral run_single \
  --outgroup Drosophila_helvetica GCA_963969585.1 \
  --species Drosophila_pseudoobscura GCF_009870125.1 \
           Drosophila_miranda GCF_003369915.1 \
  --output ./output \
  --mapq 60 \
  --suffix MAPQ60

---

### Multi-species analysis from a species list
coral run_multi \
  --species-list '[["Drosophila_melanogaster","GCF_000001215.4"],
                   ["Drosophila_sechellia","GCF_000006755.1"],
                   ["Drosophila_mauritiana","GCF_004382145.1"],
                   ["Drosophila_santomea","GCF_016746245.2"],
                   ["Drosophila_simulans","GCF_000754195.3"]]' \
  --outgroup Drosophila_simulans \
  --output ./output \
  --run-id drosophila_multi \
  --mapq 60

---

### Multi-species analysis from a Newick tree
coral run_multi \
  --newick-tree species_tree.nwk \
  --outgroup Drosophila_santomea \
  --output ./output

---

## Important Note on Multi-Species Mode

Multi-species analyses in CORAL are experimental.

While the method extends the reference-anchored framework to internal branches using phylogenetic structure, this mode has not yet undergone large-scale validation comparable to the pairwise sister-taxa pipeline.

Results from run_multi should therefore be interpreted with caution and are best used for exploratory analyses or method development.

---

## Output Summary

Each run produces a self-contained directory with:
- *.pileup.gz – multi-taxa pileup
- *_mutations.csv.gz – full mutation list
- *_mutations.json – mutation context counts
- Tables/*.tsv – normalized and collapsed spectra
- Plots/*.png – mutation spectra, genomic distributions, MAPQ diagnostics

---

## Caching and Cleanup

- Intermediate files are cached by default
- Use --no-cache to force recomputation
- Temporary FASTQs, BAMs, pileups, and interval files can be removed automatically

Manual cleanup is supported via:
cleanup_output.sh output/ --bams --pileup --intervals --genomes

---

## Scalability

- Each species-to-reference alignment is independent
- Runs can be parallelized across species or clades
- No full multiple sequence alignment required
- Suitable for large comparative datasets and cluster execution

---

## Installation and Environment

CORAL does not install or configure external bioinformatics tools.
Users must ensure required tools are available in their environment and in PATH.

### Required External Tools
- NCBI Datasets CLI  
  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

- SAMtools  
  https://www.htslib.org/doc/samtools.html

- BWA (default aligner)  
  https://github.com/bwa-mem2/bwa-mem2

Optional:
- PHYLIP  
  http://evolution.genetics.washington.edu/phylip.html

---

## Environment Setup (Recommended)

CORAL should be installed in a dedicated Python environment.
The recommended approach is conda, since it cleanly handles both Python and bioinformatics binaries.

### Option 1: Conda Environment (Recommended)

**1. Create a new environment**
```bash
conda create -n coral-env python=3.10 -y
conda activate coral-env
```

**2. Install required bioinformatics tools**
```bash
conda install -c conda-forge ncbi-datasets-cli -y
conda install -c bioconda samtools bwa -y
```

**Optional aligners:**
```bash
conda install -c bioconda minimap2 bbmap -y
```

**3. Install CORAL**

From the CORAL repository root:
```bash
pip install -e .
```

**4. Verify installation**
```bash
coral --help
python -c "import coral; print('CORAL OK')"
samtools --version
bwa
datasets --version
```

### Option 2: Python venv (Advanced users)

Use this only if external tools are already installed system-wide.

```bash
python3.10 -m venv coral-env
source coral-env/bin/activate  # On Windows: coral-env\Scripts\activate
pip install -e .
```

**Note:** You must ensure that `samtools`, `bwa`, and `datasets` are available in your PATH.

---

## PHYLIP Integration (Optional)

coral run_phylip \
  --df ./output/run_id/matching_bases.csv.gz \
  --tree ./output/run_id/annotated_tree.nwk \
  --mapping ./output/run_id/species_mapping.json \
  --phylip-exe-dir /path/to/phylip/exe

---

## Using CORAL as a Python Library

from coral import MutationExtractionPipeline

pipeline = MutationExtractionPipeline(
    species_list=[("SpeciesA", "ACC1"), ("SpeciesB", "ACC2")],
    outgroup=("Outgroup", "OUT_ACC"),
    base_output_dir="./output",
)
pipeline.run()
