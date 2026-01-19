# CORAL Output Format and Naming Conventions

This document describes the output file structure and naming conventions used by CORAL.

## Output Directory Structure

CORAL creates a self-contained output directory for each run. The directory structure follows this pattern:

```
<output_dir>/
  └── <run_id>/
      ├── *.fasta                    # Genome FASTA files
      ├── *.pileup.gz               # Multi-taxa pileup file
      ├── *_mutations.csv.gz        # Full mutation lists (one per species pair)
      ├── *_mutations.json          # Mutation context counts (one per species pair)
      ├── Mutations/                 # Mutation files directory
      ├── Triplets/                  # Trinucleotide context files
      ├── Tables/                    # Normalized spectra tables
      ├── Plots/                     # Visualization plots
      ├── Intervals/                 # Read interval files (for coverage plots)
      └── pipeline_timings.json      # Pipeline execution timing information
```

## Run ID Naming Convention

The `run_id` determines the output directory name and is used in many file names.

### For `coral run_single`:

**Default run_id:**
```
<outgroup_name>__<species1_name>__<species2_name>
```

**Example:**
- Outgroup: `Saccharomyces_mikatae_IFO_1815`
- Species 1: `Saccharomyces_paradoxus`
- Species 2: `Saccharomyces_cerevisiae_S288C`
- **Run ID:** `Saccharomyces_mikatae_IFO_1815__Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C`

**With `--suffix`:**
If you specify `--suffix test`, the run_id becomes:
```
<outgroup_name>__<species1_name>__<species2_name>_test
```

### For `coral run_multi`:

**Default run_id:**
```
multi_species_run
```

**Custom run_id:**
If you specify `--run-id <custom_id>`, that value is used directly.

## File Naming Conventions

### Pileup Files

**Pattern:**
```
<run_id>.pileup.gz
```

**Example:**
```
Saccharomyces_mikatae_IFO_1815__Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C.pileup.gz
```

### Mutation Files

Mutation files are named based on a species pair and a reference genome. **Each file contains the mutations inferred to have occurred on the phylogenetic branch leading to the first listed taxon (`<taxon1>`) since its divergence from `<taxon2>`, using `<reference>` as the reference genome.** In other words, for the file `<taxon1>__<taxon2>__<reference>__mutations.csv.gz`, the mutations listed are those on the branch leading to `<taxon1>`, relative to the common ancestor with `<taxon2>`.

**Pattern:**
```
<taxon1>__<taxon2>__<reference>__mutations.csv.gz
<taxon1>__<taxon2>__<reference>__mutations.json
```

**Example:**
For a run with:
- Reference: `Saccharomyces_mikatae_IFO_1815`
- Taxon 1: `Saccharomyces_paradoxus`
- Taxon 2: `Saccharomyces_cerevisiae_S288C`

Two mutation files are created for each direction (i.e., branch), corresponding to the mutations on the branch leading to each taxon:
```
Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz      # Mutations on the S. paradoxus branch
Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.json
Saccharomyces_cerevisiae_S288C__Saccharomyces_paradoxus__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz      # Mutations on the S. cerevisiae branch
Saccharomyces_cerevisiae_S288C__Saccharomyces_paradoxus__Saccharomyces_mikatae_IFO_1815__mutations.json
```
*For example, `Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz` contains all mutations inferred to have occurred on the branch leading to `Saccharomyces_paradoxus` after its split from `Saccharomyces_cerevisiae_S288C`, using `Saccharomyces_mikatae_IFO_1815` as the outgroup/reference.*

**Location:**
- CSV files: `<run_id>/Mutations/`
- JSON files: `<run_id>/Mutations/`

### Triplet Files

**Pattern:**
```
<taxon1>__<taxon2>__<reference>__triplets.json
```

**Location:**
`<run_id>/Triplets/`

### Normalized Spectra Tables

**Pattern:**
The following files are created in the `Tables/` directory:

- `normalized_scaled.tsv` - Normalized mutation spectra (scaled)
- `collapsed_mutations.tsv` - Raw mutation counts (collapsed)
- `scaled_raw.tsv` - Scaled raw mutation counts
- `triplets.tsv` - Trinucleotide context counts

**Location:**
`<run_id>/Tables/`

### Interval Files

**Pattern:**
```
<base_bam_name>_intervals.tsv.gz
```

Where `<base_bam_name>` is derived from the BAM file name (without `.bam` extension).

**Example:**
```
Saccharomyces_paradoxus_to_Saccharomyces_mikatae_IFO_1815_intervals.tsv.gz
```

**Location:**
`<run_id>/Intervals/`

### Plot Files

**Pattern:**
Plots are generated with descriptive names:

- Mutation spectra: `*_normalized.png`, `*_raw.png`, `*_triplets.png`
- Coverage plots: `coverage_<chromosome>.png`
- Mutation density: `mutation_density_<chromosome>.png`
- MAPQ histograms: `<species>_to_<reference>_MAPQ.png`

**Location:**
`<run_id>/Plots/`

### Multi-Species Pipeline Files

For `coral run_multi`, additional files are created:

- `matching_bases.csv.gz` - Mutation matrix for phylogenetic analysis
- `annotated_tree.nwk` - Newick tree with branch annotations
- `species_mapping.json` - Mapping between species names and internal IDs
- `mutation_spectras.tsv` - Mutation spectra summary

**Location:**
`<run_id>/`

## File Formats

### Pileup Files (`.pileup.gz`)
- Format: Gzipped text file
- Content: Multi-taxa pileup format from samtools mpileup
- Contains: Reference and aligned species base calls at each position

### Mutation CSV Files (`.csv.gz`)
- Format: Gzipped CSV
- Columns: `chromosome`, `position`, `reference_base`, `taxon1_base`, `taxon2_base`, `context`, etc.
- Contains: Full list of detected mutations with genomic positions

### Mutation JSON Files (`.json`)
- Format: JSON
- Content: Mutation context counts (trinucleotide substitution counts)
- Structure: Dictionary mapping mutation types to counts

### TSV Tables (`.tsv`)
- Format: Tab-separated values
- Content: Normalized mutation spectra, collapsed counts, or triplet frequencies
- Index: Mutation types (e.g., `A>T`, `C>G` in context)

### Interval Files (`.tsv.gz`)
- Format: Gzipped TSV
- Columns: `chromosome`, `start`, `end`
- Content: Genomic intervals covered by aligned reads

## Notes

1. **Species name format**: Species names use underscores (`_`) instead of spaces, matching NCBI naming conventions.

2. **Separator**: Double underscores (`__`) are used to separate species names in file names to avoid ambiguity.

3. **Compression**: Large files (pileups, mutation CSVs, intervals) are gzipped to save space.

4. **Directory organization**: Files are organized into subdirectories (`Mutations/`, `Tables/`, `Plots/`, etc.) for clarity.

5. **Caching**: If files already exist and `--no-cache` is not specified, CORAL will skip regeneration and use existing files.

