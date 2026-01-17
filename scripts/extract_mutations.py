#!/usr/bin/env python3
import os
import sys
import argparse
import json
import gzip
import csv
from collections import defaultdict

# Constants
REMOVE_CHARS = str.maketrans('', '', '^$[]')
CHR_IDX, POSITION_IDX, REF_NUC_IDX, N_READS_1_IDX, NUC_1_IDX, N_READS_2_IDX, NUC_2_IDX = range(7)
PREV_IDX, CUR_IDX, NEXT_IDX = 0, 1, 2
REF_IDX, TAXA1_IDX, TAXA2_IDX = 0, 1, 2

def parse_args():
    parser = argparse.ArgumentParser(description="Extract mutation triplets from a pileup file.")
    parser.add_argument("reference", help="Reference species name")
    parser.add_argument("taxa1", help="Taxa 1 name")
    parser.add_argument("taxa2", help="Taxa 2 name")
    parser.add_argument("--pileup-dir", help="Directory containing the pileup.gz file")
    parser.add_argument("--output-dir", help="Directory to write summary JSONs")
    parser.add_argument("--no-full-mutations", action="store_true", help="Don't write full gzipped mutation CSVs per species")
    parser.add_argument("--no-cache", action="store_true", help="Force regeneration of outputs")
    return parser.parse_args()

def get_nuc(nuc_field):
    cleaned = nuc_field.translate(REMOVE_CHARS)
    return cleaned[0].upper() if cleaned else 'N'

def extract_triplets(fields_list):
    context = [[], [], []]
    for fields in fields_list:
        ref_nuc = get_nuc(fields[REF_NUC_IDX])
        nuc1 = get_nuc(fields[NUC_1_IDX])
        nuc2 = get_nuc(fields[NUC_2_IDX])
        context[REF_IDX].append(ref_nuc)
        context[TAXA1_IDX].append(nuc1 if nuc1 not in {',', '.'} else ref_nuc)
        context[TAXA2_IDX].append(nuc2 if nuc2 not in {',', '.'} else ref_nuc)
    return context

def detect_mutation_triplet(triplets):
    t1_mut, t2_mut = None, None
    if triplets[REF_IDX][PREV_IDX] == triplets[TAXA1_IDX][PREV_IDX] == triplets[TAXA2_IDX][PREV_IDX] and \
       triplets[REF_IDX][NEXT_IDX] == triplets[TAXA1_IDX][NEXT_IDX] == triplets[TAXA2_IDX][NEXT_IDX]:

        ref_base = triplets[REF_IDX][CUR_IDX]
        context = triplets[REF_IDX][PREV_IDX] + ref_base + triplets[REF_IDX][NEXT_IDX]

        if ref_base == triplets[TAXA1_IDX][CUR_IDX] and ref_base != triplets[TAXA2_IDX][CUR_IDX]:
            t2_mut = f"{context[0]}[{ref_base}>{triplets[TAXA2_IDX][CUR_IDX]}]{context[2]}"
        elif ref_base == triplets[TAXA2_IDX][CUR_IDX] and ref_base != triplets[TAXA1_IDX][CUR_IDX]:
            t1_mut = f"{context[0]}[{ref_base}>{triplets[TAXA1_IDX][CUR_IDX]}]{context[2]}"
    return t1_mut, t2_mut

def quality_check(fields):
    return fields and '*' not in fields[NUC_1_IDX] and '*' not in fields[NUC_2_IDX] and \
            all_same(fields[NUC_1_IDX].translate(REMOVE_CHARS)) and all_same(fields[NUC_2_IDX].translate(REMOVE_CHARS))

def all_same(seq):
    return len(seq) > 0 and all(ch == seq[0] for ch in seq)

def parse_line(line):
    parts = line.strip().split('\t')
    return parts[:5] + parts[6:-1] if len(parts) >= 9 else None

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pileup_file = os.path.join(args.pileup_dir, f"{args.reference}__{args.taxa1}__{args.taxa2}.pileup.gz")
    out_json1 = os.path.join(args.output_dir, f"{args.taxa1}__{args.taxa2}__{args.reference}__mutations.json")
    out_json2 = os.path.join(args.output_dir, f"{args.taxa2}__{args.taxa1}__{args.reference}__mutations.json")

    csv_path1 = os.path.join(args.output_dir, f"{args.taxa1}__{args.taxa2}__{args.reference}__mutations.csv.gz") if not args.no_full_mutations else None
    csv_path2 = os.path.join(args.output_dir, f"{args.taxa2}__{args.taxa1}__{args.reference}__mutations.csv.gz") if not args.no_full_mutations else None
    
    jsons_exist = all(os.path.exists(p) for p in [out_json1, out_json2])
    csvs_exist = not args.no_full_mutations and all(os.path.exists(p) for p in [csv_path1, csv_path2])

    if not args.no_cache and jsons_exist and (args.no_full_mutations or csvs_exist):
        print("Mutation counts already exist. Skipping.")
        return

    species_mut1 = defaultdict(int)
    species_mut2 = defaultdict(int)

    writer1 = gzip.open(csv_path1, 'wt', newline='') if csv_path1 else None
    writer2 = gzip.open(csv_path2, 'wt', newline='') if csv_path2 else None
    csv_writer1 = csv.writer(writer1) if writer1 else None
    csv_writer2 = csv.writer(writer2) if writer2 else None

    if csv_writer1:
        csv_writer1.writerow(["chromosome", "position", "reference_base", "mutated_base", "context", "triplet"])
    if csv_writer2:
        csv_writer2.writerow(["chromosome", "position", "reference_base", "mutated_base", "context", "triplet"])

    with gzip.open(pileup_file, 'rt') as f:
        line_fields = [None, parse_line(f.readline()), parse_line(f.readline())]
        qc_flags = [False, quality_check(line_fields[1]), quality_check(line_fields[2])]

        for line in f:
            line_fields = [line_fields[1], line_fields[2], parse_line(line)]
            qc_flags = [qc_flags[1], qc_flags[2], quality_check(line_fields[2])]

            if all(qc_flags):
                triplets = extract_triplets(line_fields)
                mut1, mut2 = detect_mutation_triplet(triplets)
                chrom = line_fields[1][CHR_IDX]
                pos = int(line_fields[1][POSITION_IDX])
                context = ''.join(triplets[REF_IDX])
                ref_base = triplets[REF_IDX][CUR_IDX]

                if mut1:
                    species_mut1[mut1] += 1
                    if csv_writer1:
                        csv_writer1.writerow([chrom, pos, ref_base, triplets[TAXA1_IDX][CUR_IDX], context, mut1])
                if mut2:
                    species_mut2[mut2] += 1
                    if csv_writer2:
                        csv_writer2.writerow([chrom, pos, ref_base, triplets[TAXA2_IDX][CUR_IDX], context, mut2])

    if writer1:
        writer1.close()
        print(f"Full mutation CSV written to: {csv_path1}")
    if writer2:
        writer2.close()
        print(f"Full mutation CSV written to: {csv_path2}")

    with open(out_json1, 'w') as f:
        json.dump(species_mut1, f, indent=2)
    with open(out_json2, 'w') as f:
        json.dump(species_mut2, f, indent=2)
    print(f"Saved mutation counts to {out_json1} and {out_json2}")

if __name__ == "__main__":
    main()
