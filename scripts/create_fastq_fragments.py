#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO
from pathlib import Path

def fastq_fragments(input_fasta, output_fastq, length=150, offset=75, force=False):
    if os.path.exists(output_fastq) and not force:
        print(f"{output_fastq} already exists. Skipping.")
        return

    chromosomes = {
        record.id: record.seq
        for record in SeqIO.parse(input_fasta, "fasta")
    }

    def split_sequence(seq, offset, section_size):
        last_start = len(seq) - section_size
        sections = [(i, seq[i:i+section_size]) for i in range(0, last_start, offset)]
        if last_start > 0:
            sections.append((last_start, seq[last_start:]))
        return sections

    with open(output_fastq, 'w') as out:
        for chrom_name, sequence in chromosomes.items():
            for i, (start, frag) in enumerate(split_sequence(sequence, offset, length)):
                end = start + len(frag) - 1
                out.write(f"@{chrom_name}_{start+1}_{end+1}\n")
                out.write(f"{frag}\n+\n")
                out.write(f"{'I' * len(frag)}\n")

    print(f"Wrote {output_fastq}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create synthetic FASTQ reads from a genome FASTA.")
    parser.add_argument("input_fasta", help="Path to the input genome FASTA file")
    parser.add_argument("--output", "-o", help="Path to output FASTQ file (default: same dir, .fastq extension)")
    parser.add_argument("--length", type=int, default=150, help="Length of each fragment (default: 150)")
    parser.add_argument("--offset", type=int, default=75, help="Step size between fragments (default: 75)")
    parser.add_argument("--force", action="store_true", help="Overwrite output file if it exists")
    args = parser.parse_args()

    input_path = Path(args.input_fasta)

    # Determine output path
    if args.output:
        output_fastq = Path(args.output)
    else:
        output_fastq = input_path.with_suffix(".fastq")

    fastq_fragments(
        input_fasta=str(input_path),
        output_fastq=str(output_fastq),
        length=args.length,
        offset=args.offset,
        force=args.force
    )
