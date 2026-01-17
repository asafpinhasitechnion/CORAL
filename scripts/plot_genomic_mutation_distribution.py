import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


def parse_fai(fai_file):
    chrom_lengths = {}
    with open(fai_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            length = int(fields[1])
            chrom_lengths[chrom] = length
    return chrom_lengths


def compute_mutation_density(mutation_file, chrom, chrom_length, bin_size=1000, slide=1000, mut_regex=None):
    df = pd.read_csv(mutation_file, compression='infer')
    df = df[df['chromosome'] == chrom]
    if mut_regex:
        df = df[df["triplet"].str.contains(mut_regex, regex=True, na=False)]

    starts = np.arange(0, chrom_length - bin_size + 1, slide)
    mutation_counts = []

    for bin_start in starts:
        bin_end = bin_start + bin_size
        count = ((df['position'] >= bin_start) & (df['position'] < bin_end)).sum()
        mutation_counts.append(count)

    midpoints = [start + bin_size // 2 for start in starts]
    return midpoints, mutation_counts


def compute_binned_coverage(interval_file, chrom, chrom_length, bin_size=1000, slide=1000, mut_regex=None):
    df = pd.read_csv(
        interval_file,
        sep='\t',
        compression='infer',
        dtype={"chromosome": str, "start": int, "end": int},
        header=0  # assumes the file includes a header line
    )
    df[df["chromosome"] == chrom].sort_values("start").reset_index(drop=True)

    starts = np.arange(0, chrom_length - bin_size + 1, slide)
    coverage = []

    current_index = 0
    n = len(df)

    for bin_start in starts:
        bin_end = bin_start + bin_size
        total_overlap = 0

        while current_index < n and df.at[current_index, "end"] <= bin_start:
            current_index += 1

        check_index = current_index
        while check_index < n and df.at[check_index, "start"] < bin_end:
            read_start = df.at[check_index, "start"]
            read_end = df.at[check_index, "end"]
            overlap_start = max(read_start, bin_start)
            overlap_end = min(read_end, bin_end)
            total_overlap += max(0, overlap_end - overlap_start)
            check_index += 1

        coverage.append(total_overlap / bin_size)

    return coverage


def plot_mutation_density(midpoints_list, mutation_list, labels, chrom, output_path, normalized=False, mut_regex=None):
    plt.figure(figsize=(15, 5))
    for midpoints, values, label in zip(midpoints_list, mutation_list, labels):
        plt.plot(midpoints, values, label=label, lw=1.5)

    plt.title(f"Mutation Density over {chrom}{' (normalized)' if normalized else ''}")
    plt.xlabel("Genomic Position")
    plt.ylabel("Mutation Density")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot mutation density across a chromosome, optionally normalized by coverage.")
    parser.add_argument("--mutation-dir", required=True, help="Directory with mutation .csv.gz files")
    parser.add_argument("--fai-file", required=True, help=".fai index file for chromosome lengths")
    parser.add_argument("--chromosome", required=True, help="Chromosome to plot")
    parser.add_argument("--output-dir", required=True, help="Output directory path")
    parser.add_argument("--bin-size", type=int, default=100000)
    parser.add_argument("--slide", type=int)
    parser.add_argument("--coverage-dir", help="Optional directory of interval files for normalization")
    parser.add_argument("--mutation-category", help="Regex for a mutation type to include")
    args = parser.parse_args()

    if args.slide is None:
        args.slide = args.bin_size

    chrom_lengths = parse_fai(args.fai_file)
    chrom = args.chromosome
    chrom_length = chrom_lengths.get(chrom)
    if chrom_length is None:
        raise ValueError(f"Chromosome {chrom} not found in {args.fai_file}")

    mutation_files = [os.path.join(args.mutation_dir, f) for f in os.listdir(args.mutation_dir) if f.endswith(".csv.gz")]
    mutation_files.sort()
    labels = [os.path.basename(f).replace("_mutations.csv.gz", "") for f in mutation_files]

    if args.coverage_dir:
        interval_files = [
            os.path.join(args.coverage_dir, f)
            for f in os.listdir(args.coverage_dir)
            if f.endswith(".tsv") or f.endswith(".tsv.gz")
        ]
        interval_files.sort()
    else:
        interval_files = []

    mutation_regex = re.compile(args.mutation_category) if args.mutation_category else None

    midpoints_list = []
    values_list = []

    for i, mutation_file in enumerate(mutation_files):
        midpoints, mutations = compute_mutation_density(
            mutation_file, chrom, chrom_length, args.bin_size, args.slide, mutation_regex
        )

        if args.coverage_dir:
            coverage_file = interval_files[i]
            coverage = compute_binned_coverage(
                coverage_file, chrom, chrom_length, args.bin_size, args.slide, mutation_regex
            )
            mutations = [m / c if c > 0 else 0 for m, c in zip(mutations, coverage)]

        midpoints_list.append(midpoints)
        values_list.append(mutations)

    # Create meaningful plot name from arguments
    plot_name = f"mutation_density_{chrom}"
    if args.mutation_category:
        plot_name += f"_{args.mutation_category}"
    if args.coverage_dir:
        plot_name += "_normalized"
    plot_name += ".png"

    output_path = os.path.join(args.output_dir, plot_name)
    plot_mutation_density(midpoints_list, values_list, labels, chrom, output_path, normalized=bool(args.coverage_dir), mut_regex=mutation_regex)


if __name__ == "__main__":
    main()
