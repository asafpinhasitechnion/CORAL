#!/usr/bin/env python3
import os
import sys
import argparse
import json
import gzip
from collections import defaultdict

REMOVE_CHARS = str.maketrans('', '', '^$[]')
CHR_IDX, POSITION_IDX, REF_NUC_IDX, N_READS_1_IDX, NUC_1_IDX, N_READS_2_IDX, NUC_2_IDX = range(7)
PREV_IDX, CUR_IDX, NEXT_IDX = 0, 1, 2
REF_IDX, TAXA1_IDX, TAXA2_IDX = 0, 1, 2

def parse_args():
    parser = argparse.ArgumentParser(description="Extract 3-mer triplet frequencies from a pileup file.")
    parser.add_argument("reference", help="Reference species name")
    parser.add_argument("taxa1", help="Taxa 1 name")
    parser.add_argument("taxa2", help="Taxa 2 name")
    parser.add_argument("--pileup-dir", help="Input directory containing .pileup.gz")
    parser.add_argument("--output-dir", help="Output directory for JSONs")
    parser.add_argument("--no-cache", action="store_true", help="Force recomputation even if outputs exist")
    return parser.parse_args()

def all_same(seq):
    return len(seq) > 0 and all(ch == seq[0] for ch in seq)

def get_nuc(field):
    cleaned = field.translate(REMOVE_CHARS)
    return cleaned[0].upper() if cleaned else 'N'


def parse_line(line):
    parts = line.strip().split('\t')
    return parts[:5] + parts[6:-1] if len(parts) >= 9 else None

def passes_qc(fields):
    return fields is not None and \
           '*' not in fields[NUC_1_IDX] and '*' not in fields[NUC_2_IDX] and \
           all_same(fields[NUC_1_IDX].translate(REMOVE_CHARS)) and \
           all_same(fields[NUC_2_IDX].translate(REMOVE_CHARS))

def extract_triplets(line_fields):
    sequences = [[], [], []]
    for fields in line_fields:
        ref_nuc = get_nuc(fields[REF_NUC_IDX])
        nuc1 = get_nuc(fields[NUC_1_IDX])
        nuc2 = get_nuc(fields[NUC_2_IDX])
        nuc1 = ref_nuc if nuc1 in {',', '.'} else nuc1
        nuc2 = ref_nuc if nuc2 in {',', '.'} else nuc2
        sequences[REF_IDX].append(ref_nuc)
        sequences[TAXA1_IDX].append(nuc1)
        sequences[TAXA2_IDX].append(nuc2)
    return sequences

def relevant_triplet(triplets):
    if triplets[REF_IDX][PREV_IDX] == triplets[TAXA1_IDX][PREV_IDX] == triplets[TAXA2_IDX][PREV_IDX] and \
       triplets[REF_IDX][NEXT_IDX] == triplets[TAXA1_IDX][NEXT_IDX] == triplets[TAXA2_IDX][NEXT_IDX]:
        if triplets[REF_IDX][CUR_IDX] == triplets[TAXA1_IDX][CUR_IDX] or \
            triplets[REF_IDX][CUR_IDX] == triplets[TAXA2_IDX][CUR_IDX]:
            return True
    return False

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pileup_file = os.path.join(args.pileup_dir, f"{args.reference}__{args.taxa1}__{args.taxa2}.pileup.gz")
    out_json1 = os.path.join(args.output_dir, f"{args.taxa1}__{args.taxa2}__{args.reference}__triplets.json")
    out_json2 = os.path.join(args.output_dir, f"{args.taxa2}__{args.taxa1}__{args.reference}__triplets.json")

    if all(os.path.exists(p) for p in [out_json1, out_json2]) and not args.no_cache:
        print("Triplet counts already exist. Skipping.")
        return

    triplet_dict1 = defaultdict(int)
    triplet_dict2 = defaultdict(int)

    with gzip.open(pileup_file, 'rt') as f:
        line_fields = [None, parse_line(f.readline()), parse_line(f.readline())]
        qc_flags = [False, passes_qc(line_fields[1]), passes_qc(line_fields[2])]

        for line in f:
            line_fields = [line_fields[CUR_IDX], line_fields[NEXT_IDX], parse_line(line)]
            qc_flags = [qc_flags[CUR_IDX], qc_flags[NEXT_IDX], passes_qc(line_fields[NEXT_IDX])]
            if all(qc_flags):
                triplets = extract_triplets(line_fields)
                # triplet_dict1[''.join(triplets[TAXA1_IDX])] += 1
                # triplet_dict2[''.join(triplets[TAXA2_IDX])] += 1 
                if relevant_triplet(triplets):
                    triplet = ''.join(triplets[REF_IDX])
                    triplet_dict1[triplet] += 1
                    triplet_dict2[triplet] += 1


    with open(out_json1, 'w') as f:
        json.dump(triplet_dict1, f, indent=2)
    with open(out_json2, 'w') as f:
        json.dump(triplet_dict2, f, indent=2)

    print(f"Triplet dictionaries written to {out_json1} and {out_json2}")

if __name__ == "__main__":
    main()
