import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
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


def compute_binned_coverage(interval_file, chrom, chrom_length, bin_size=1000, slide=1000):
    # Support .tsv and .tsv.gz
    df = pd.read_csv(
        interval_file,
        sep='\t',
        compression='infer',
        dtype={"chrom": str, "start": int, "end": int},
        header=0  # assumes the file includes a header line
    )

    df = df[df["chromosome"] == chrom].sort_values("start").reset_index(drop=True)

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

    midpoints = [start + bin_size // 2 for start in starts]
    return midpoints, coverage


def plot_coverage(midpoints_list, coverage_list, labels, chrom, output_path):
    plt.figure(figsize=(15, 5))
    for midpoints, coverage, label in zip(midpoints_list, coverage_list, labels):
        plt.plot(midpoints, coverage, label=label, lw=1.5)

    plt.title(f"Coverage over {chrom} — {' vs '.join(labels)}")
    plt.xlabel("Genomic Position")
    plt.ylabel("Normalized Read Coverage (per base)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Generate coverage plot for a given chromosome using interval files.")
    parser.add_argument("--interval-dir", required=True, help="Directory containing interval .tsv[.gz] files")
    parser.add_argument("--fai-file", required=True, help="FASTA index file (.fai) with chromosome lengths")
    parser.add_argument("--chromosome", required=True, help="Chromosome to plot")
    parser.add_argument("--output", required=True, help="Path to save the output plot")
    parser.add_argument("--bin-size", type=int, default=100000, help="Size of each bin")
    parser.add_argument("--slide", type=int, help="Sliding window size")
    args = parser.parse_args()

    if args.slide is None:
        args.slide = args.bin_size

    chrom_lengths = parse_fai(args.fai_file)
    chrom = args.chromosome
    chrom_length = chrom_lengths.get(chrom)

    if chrom_length is None:
        raise ValueError(f"Chromosome {chrom} not found in {args.fai_file}")

    # Support both .tsv and .tsv.gz files
    interval_files = [
        os.path.join(args.interval_dir, f)
        for f in os.listdir(args.interval_dir)
        
        if f.endswith(".tsv") or f.endswith(".tsv.gz")
    ]
    labels = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in interval_files]

    midpoints_list = []
    coverage_list = []

    for file in interval_files:
        midpoints, coverage = compute_binned_coverage(file, chrom, chrom_length, args.bin_size, args.slide)
        midpoints_list.append(midpoints)
        coverage_list.append(coverage)

    plot_coverage(midpoints_list, coverage_list, labels, chrom, args.output)


if __name__ == "__main__":
    main()
