
# Multi-Species Mutation Extraction Pipeline

This repository contains a tool for scalable extraction, detection and analysis of mutations in species' evolutionary history.
The pipeline aligns multiple species to a reference genome, simulates reads, filters by mapping quality, extracts unambiguous triplet mutations, and visualizes mutation and coverage patterns.

---

## Core Concepts

### Mutation Context Extraction
The pipeline identifies positions in the genome where:
- All species share the **same flanking bases**
- One species has a **distinct middle base**

This allows the extraction of **unambiguous point mutations** in triplet format, enabling detailed mutation spectrum analysis. 
The distribution of mutation spectra, mutation rates and mutation types across the reference genome can be visualized and processed.

### Distributed Alignment Strategy Using "Pseudo-Reads"
Each non-outgroup species is aligned independently to a shared **reference genome (outgroup)**. This avoids full multiple alignment and makes the process **scalable**, while minimizing the notorious MSA reference-bias. 
The tool uses state-of-the-art local aligners, normaly used aligning DNA sequencing reads to a reference genome, optimized for efficiency.
The tool is increadibly lightweight compared to standard whole-genome multiple sequence aligners, making it **increadibly fit for parallel running** across different species groups.

---

## Functional Workflow

### Step 1: Genome Preparation
- Download genomes by NCBI accession
- Index reference genome for alignment tools

### Step 2: Read Simulation and Alignment
- Simulate FASTQ reads by sliding a window across the genome
- Align simulated reads to the outgroup
- Filter alignments by MAPQ and coverage
- Allows user to customize the aligner and its parameters

### Step 3: Mutation Detection
- Generate a pileup of reference + aligned BAMs
- Extract unambiguous triplet mutation events
- Optionally: extract the full mutations list, including genomic position

### Step 4: Normalization & Analysis
- Normalize mutation counts by underlying triplet abundance
- Collapse complementary strands for canonical spectra
- Visualize results: mutation spectra, genomic distribution, coverage

---

## Example Usage

### Run the full pipeline
```bash
python run_pipeline.py Pieris_rapae GCF_905147795.1 \
                       Pieris_napi GCF_905475465.1 \
                       Leptophobia_aripa GCA_951799465.1 \
                       --mapq 60 --aligner bwa --remove-temp --genomic-position-plots
```

### Use an individual tool
```bash
python extract_mutations.py Leptophobia_aripa Pieris_rapae Pieris_napi \
  --pileup-dir Output --output-dir Output
```

---

## Caching & Cleanup

All stages support optional caching:
- Use `--no-cache` to **force regeneration** of intermediate files
- Default behavior: **skip re-running** steps if outputs exist

To manage disk usage, pass:
```bash
--remove-temp
```
This deletes intermediate BAMs, pileups, FASTQs, and intervals once they're no longer needed.

A dedicated script `cleanup_output.sh` allows **manual cleanup** by component:
```bash
bash cleanup_output.sh Output/ --bams --pileup --intervals --genomes
```

---

## Scalability & Parallelization

- Each **species vs reference** alignment is independent → can be parallelized
- Output's for each run are saved in separate folders, allowing parallelization without conflict
- No full multiple sequence alignment required
- Easily scaled to hundreds of genomes with a job scheduler (e.g. SLURM/Condor)

This makes the pipeline ideal for large-scale evolutionary or mutational studies across clades.

---

## Output Summary

- `*.pileup.gz` – multi-taxa pileup for reference + aligned BAMs
- `*_mutations.csv.gz` – full triplet mutation calls
- `*_mutations.json` – mutation context counts
- `Tables/*.tsv` – normalized, collapsed spectra
- `Plots/*.png` – spectra, genome-wide mutation distributions, MAPQ histograms

---

## Optional Visualizations

- **Mutation spectra** (raw & normalized)
- **Triplet distributions** per taxon
- **MAPQ score histogram** (for filtering diagnostics)
- **Mutation & coverage distribution** across chromosomes

---

## Requirements & Installation

### ⚠️ **Platform Requirements**

**Linux or WSL2 (Windows Subsystem for Linux) is REQUIRED.**

This pipeline relies on Linux-specific bioinformatics tools (`samtools`, `bwa`, etc.) that are not available on native Windows. If you are using Windows:

1. **Install WSL2** (recommended):
   ```powershell
   # In Windows PowerShell (as Administrator)
   wsl --install
   ```
   Then restart your computer and use the Ubuntu terminal.

2. **Alternative**: Use a Linux virtual machine (VirtualBox, VMware, etc.)

All following installation instructions assume a Linux/WSL environment.

---

### **Required External Tools**

These tools must be installed and available in your `PATH` before running the pipeline:

#### 1. **NCBI Datasets CLI** (Required for genome downloads)

The pipeline uses the NCBI Datasets command-line tool to download genomes.

**Installation:**

```bash
# Download and install the datasets CLI
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
chmod +x datasets
sudo mv datasets /usr/local/bin/

# Verify installation
datasets --version
```

**Alternative (using conda):**
```bash
conda install -c conda-forge ncbi-datasets-cli
```

**Documentation:** https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

---

#### 2. **SAMtools** (Required for BAM processing and pileup generation)

**Installation using conda (recommended):**
```bash
conda install -c bioconda samtools=1.17
```

**Or using package manager (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install samtools
```

**Or build from source:**
```bash
# Install dependencies
sudo apt-get install libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses-dev

# Download and compile
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
tar -xjf samtools-1.17.tar.bz2
cd samtools-1.17
./configure --prefix=/usr/local
make
sudo make install
```

**Verify installation:**
```bash
samtools --version
```

---

#### 3. **BWA** (Required aligner - default)

**Installation using conda (recommended):**
```bash
conda install -c bioconda bwa=0.7.17
```

**Or using package manager (Ubuntu/Debian):**
```bash
sudo apt-get install bwa
```

**Or build from source:**
```bash
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
sudo mv bwa /usr/local/bin/
```

**Verify installation:**
```bash
bwa
```

---

#### 4. **Standard Unix Tools** (Usually pre-installed)

- `gzip` - For compression (usually pre-installed)
- `unzip` - For extracting downloaded genome archives

**If missing, install:**
```bash
sudo apt-get install gzip unzip
```

---

#### 5. **PHYLIP** (Optional - Required for `run_phylip` command)

PHYLIP is needed for phylogenetic analysis of mutation matrices. The pipeline expects PHYLIP executables (e.g., `dnapars`, `dnapenny`) in a specific directory.

**Option A: Install via package manager (Easiest)**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install phylip

# The executables will be in /usr/bin/phylip/ or similar
# You'll need to specify this path using --phylip-exe-dir
```

**Option B: Download and compile PHYLIP (Recommended for control)**

1. **Download PHYLIP:**
   ```bash
   cd /tmp
   wget https://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz
   tar -xzf phylip-3.697.tar.gz
   cd phylip-3.697/src
   ```

2. **Compile:**
   ```bash
   make -f Makefile.unx install
   # This will compile all PHYLIP programs
   ```

3. **Locate executables:**
   The compiled executables will be in `phylip-3.697/exe/`. Common programs include:
   - `dnapars` - DNA parsimony
   - `dnapenny` - DNA branch-and-bound
   - `dnamlk` - DNA maximum likelihood

4. **Set up for CORAL:**
   
   **Option 1:** Place PHYLIP executables in the default location expected by CORAL:
   ```bash
   # Copy PHYLIP directory to CORAL source
   cp -r phylip-3.697 /path/to/CORAL/src/
   ```

   **Option 2:** Use `--phylip-exe-dir` when running PHYLIP commands:
   ```bash
   coral run_phylip --phylip-exe-dir /path/to/phylip-3.697/exe ...
   ```

**Verify installation:**
```bash
# Test if dnapars is available
/path/to/phylip-3.697/exe/dnapars
# (Press Ctrl+C after it starts - this just verifies it runs)
```

**PHYLIP Documentation:** http://evolution.genetics.washington.edu/phylip.html

---

### **Optional Aligners**

These can be used instead of BWA by specifying `--aligner-name`:

- **minimap2**: `conda install -c bioconda minimap2`
- **bbmap**: `conda install -c bioconda bbmap`

---

### **Python Dependencies**

#### Installation Options:

**Option 1: Using conda (Recommended - includes bioinformatics tools)**
```bash
# Create environment from file
conda env create -f environment.yml
conda activate coral-env

# Install CORAL package
pip install -e .
```

**Option 2: Using pip only**
```bash
# Create virtual environment
python3.10 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install package (dependencies will be installed automatically)
pip install -e .
```

**Python Dependencies (automatically installed):**
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `matplotlib` - Plotting
- `biopython` - Bioinformatics utilities
- `pysam` - BAM/SAM file handling
- `ete3` - Phylogenetic tree handling
- `psutil` - System and process utilities
- `click` - CLI framework

**Verify Python installation:**
```bash
python --version  # Should be ≥ 3.10
pip list  # Check all dependencies are installed
```

---

### **Installation Verification Checklist**

Before running the pipeline, verify all tools are installed:

```bash
# Check external tools
datasets --version      # NCBI Datasets CLI
samtools --version      # SAMtools
bwa                     # BWA (should print help)
gzip --version          # gzip
unzip -v                # unzip

# Check Python
python --version        # Should be ≥ 3.10

# Check Python package
python -c "import coral; print('CORAL installed successfully')"

# Check PHYLIP (if using run_phylip)
ls /path/to/phylip-3.697/exe/dnapars  # Should exist
```

---

### **Quick Start Installation (All-in-one)**

For a complete setup on Ubuntu/WSL, you can run:

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y build-essential wget curl unzip gzip

# Install NCBI Datasets CLI
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
chmod +x datasets
sudo mv datasets /usr/local/bin/

# Install via conda (includes samtools, bwa, and Python packages)
conda env create -f environment.yml
conda activate coral-env
pip install -e .

# (Optional) Install PHYLIP
sudo apt-get install -y phylip
# OR download and compile as described above
```

---

## Usage Examples

### Single pipeline (3-species: outgroup + 2 ingroup)
```bash
coral run_single \
  --outgroup Drosophila_helvetica GCA_963969585.1 \
  --species Drosophila_pseudoobscura GCF_009870125.1 Drosophila_miranda GCF_003369915.1 \
  --output ./Output_OO \
  --mapq 60 \
  --suffix MAPQ60
```

### Multi-species pipeline from Newick tree
```bash
coral run_multi \
  --newick-tree "(((Drosophila_sechellia|GCF_004382195.2,Drosophila_melanogaster|GCF_000001215.4),Drosophila_mauritiana|GCF_004382145.1),Drosophila_santomea|GCF_016746245.2);" \
  --output ./Output_OO \
  --outgroup Drosophila_santomea
```

### Multi-species pipeline from species list
```bash
coral run_multi \
  --species-list '[["Drosophila_sechellia","GCF_004382195.2"],["Drosophila_melanogaster","GCF_000001215.4"],["Drosophila_mauritiana","GCF_004382145.1"],["Drosophila_santomea","GCF_016746245.2"]]' \
  --output ./Output_OO \
  --outgroup Drosophila_santomea \
  --run-id drosophila1_run_mutiple_species
```

### Generate plots from existing output
```bash
coral plot --tables ./output/run_id/Tables
```

### Run PHYLIP on mutation matrix
```bash
coral run_phylip \
  --df ./output/run_id/matching_bases.csv.gz \
  --tree ./output/run_id/annotated_tree.nwk \
  --mapping ./output/run_id/species_mapping.json
```

### Using as a Python library
```python
from coral import MutationExtractionPipeline

pipeline = MutationExtractionPipeline(
    species_list=[("Species1", "ACC1"), ("Species2", "ACC2")],
    outgroup=("Outgroup", "OUT_ACC"),
    base_output_dir="./output"
)
pipeline.run()
```
