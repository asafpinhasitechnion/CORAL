# CORAL:  Comparative Orthologous Read-based Analysis of Lineage Substitutions

CORAL is a tool for scalable extraction, detection, and analysis of point mutations across species evolutionary history.
It aligns multiple species to a shared reference genome, simulates reads, filters alignments by mapping quality, extracts unambiguous trinucleotide substitutions, and summarizes mutation rates and spectra.

Reference  
Please cite: will be available upon publication

---
<img width="4066" height="1176" alt="CORAL_pipeline" src="https://github.com/user-attachments/assets/dd9d9d43-8775-4585-9be7-1f0bafebfc92" />

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

Standard short-read aligners (default: BWA, validated with classic BWA; BWA-MEM2 and other aligners are also supported) are used for efficiency and robustness.

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

### Three-taxon pipeline (pairwise sister-taxa analysis: outgroup + two ingroups)
coral run_single \
  --outgroup Saccharomyces_mikatae_IFO_1815 GCF_947241705.1 \
  --species Saccharomyces_paradoxus GCF_002079055.1 Saccharomyces_cerevisiae_S288C GCF_000146045.2 \
  --output ../test_output \
  --mapq 60 \
  --suffix test

---

### Multi-species analysis from a species list
coral run_multi \
  --species-list '[["Drosophila_melanogaster","GCF_000001215.4"],["Drosophila_sechellia","GCF_000006755.1"],["Drosophila_mauritiana","GCF_004382145.1"],["Drosophila_santomea","GCF_016746245.2"],["Drosophila_simulans","GCF_016746395.2"]]' \
  --outgroup Drosophila_simulans \
  --output ../test_output \
  --run-id multi_test \
  --mapq 60

---

### Multi-species analysis from a Newick tree
coral run_multi \
  --newick-tree "(((Drosophila_melanogaster|GCF_000001215.4,Drosophila_sechellia|GCF_000006755.1),Drosophila_mauritiana|GCF_004382145.1),Drosophila_santomea|GCF_016746245.2,Drosophila_simulans|GCF_016746395.2);" \
  --outgroup Drosophila_simulans \
  --output ../test_output \
  --run-id multi_test \
  --mapq 60

---

## Important Note on Multi-Species Mode

Multi-species analyses in CORAL are experimental.

While the method extends the reference-anchored framework to internal branches using phylogenetic structure, this mode has not yet undergone large-scale validation comparable to the pairwise sister-taxa pipeline.

Results from run_multi should therefore be interpreted with caution and are best used for exploratory analyses or method development.

---

## Output Summary

Each run produces a self-contained directory with:
- `<run_id>.pileup.gz` – multi-taxa pileup (in run root)
- `Mutations/*_mutations.csv.gz` – full mutation lists (one per species pair)
- `Mutations/*_mutations.json` – mutation context counts (one per species pair)
- `Tables/*.tsv` – normalized and collapsed spectra
- `Plots/*.png` – mutation spectra, genomic distributions, MAPQ diagnostics

**Mutation file naming:** Files are named `<taxon1>__<taxon2>__<reference>__mutations.*` where the file contains mutations inferred to have occurred on the branch leading to `<taxon1>` since its divergence from `<taxon2>`, using `<reference>` as the outgroup/reference genome. See [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) for detailed naming conventions and file formats.

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

### Option 1: Using environment.yml (Recommended - Easiest)

The `environment.yml` file automatically installs all required dependencies, including external bioinformatics tools.

**1. Create environment from environment.yml**
```bash
conda env create -f environment.yml
conda activate coral-env
```

**2. Install CORAL**
```bash
pip install -e .
```

**3. Verify installation**
```bash
coral --help
python -c "import coral; print('CORAL OK')"
samtools --version
bwa
datasets --version
```

**Note:** The `environment.yml` includes:
- Python 3.10
- BWA (classic, version 0.7.17)
- SAMtools (version 1.17)
- NCBI Datasets CLI
- All Python dependencies

**For multi-species pipeline:** PHYLIP is commented out in `environment.yml`. Uncomment the `phylip` line in the file, or install separately:
```bash
conda install -c bioconda phylip -y
```

---

### Option 2: Manual Conda Installation

If you prefer to set up the environment manually:

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

**Optional: PHYLIP (Required for multi-species pipeline)**
```bash
conda install -c bioconda phylip -y
```
PHYLIP must be available in PATH for the multi-species pipeline (`coral run_multi`). When installed via conda, PHYLIP executables are automatically detected. The pipeline will raise an error if PHYLIP is not found in PATH.

**3. Install CORAL**
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
# If using multi-species pipeline, verify PHYLIP:
dnapars  # Should show PHYLIP help if installed
```

---

### Required External Tools

If installing manually (not using `environment.yml`), you need:

- **NCBI Datasets CLI**  
  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

- **SAMtools**  
  https://www.htslib.org/doc/samtools.html

- **BWA (default aligner)** - Classic BWA is the default and was used for validation. Other versions including BWA-MEM2 are also supported.  
  Classic BWA (default): https://github.com/lh3/bwa  
  BWA-MEM2: https://github.com/bwa-mem2/bwa-mem2  
  Note: `conda install bwa` installs classic BWA. To use BWA-MEM2, use `--aligner-name bwa-mem2` (BWA-MEM2 must be installed separately).

**Optional (Required for multi-species pipeline):**
- **PHYLIP** - Required for `coral run_multi` command  
  Install via conda: `conda install -c bioconda phylip`  
  http://evolution.genetics.washington.edu/phylip.html

### Option 3: Python venv (Advanced users)

Use this only if external tools are already installed system-wide.

```bash
python3.10 -m venv coral-env
source coral-env/bin/activate  # On Windows: coral-env\Scripts\activate
pip install -e .
```

**Note:** You must ensure that `samtools`, `bwa`, and `datasets` are available in your PATH.

---

## PHYLIP Integration

PHYLIP is **required for the full multi-species pipeline** (`coral run_multi`) and optional for standalone phylogenetic analysis (`coral run_phylip`).

### Installation

**Recommended (via conda):**
```bash
conda install -c bioconda phylip
```

When installed via conda, PHYLIP executables are automatically detected in PATH. No additional configuration is needed.

### Usage

**Multi-species pipeline (automatically uses PHYLIP):**
```bash
coral run_multi --newick-tree "..." --outgroup SpeciesName --output ./output
```

**Standalone PHYLIP analysis:**
```bash
coral run_phylip \
  --df ./output/run_id/matching_bases.csv.gz \
  --tree ./output/run_id/annotated_tree.nwk \
  --mapping ./output/run_id/species_mapping.json
```

**Note:** PHYLIP must be installed and available in PATH. The multi-species pipeline will raise an error if PHYLIP is not found. Conda installation is recommended but not required.

---

## Documentation

- **[tutorial.ipynb](tutorial.ipynb)** - Command-line tutorial with examples
- **[OUTPUT_FORMAT.md](OUTPUT_FORMAT.md)** - Output file structure and naming conventions

---

## Using CORAL as a Python Library

from coral import MutationExtractionPipeline

pipeline = MutationExtractionPipeline(
    species_list=[("SpeciesA", "ACC1"), ("SpeciesB", "ACC2")],
    outgroup=("Outgroup", "OUT_ACC"),
    base_output_dir="./output",
)
pipeline.run()
