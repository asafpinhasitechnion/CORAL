#!/usr/bin/env python3
import os
import argparse
import json
import gzip
from collections import defaultdict

# Constants
REMOVE_CHARS = str.maketrans('', '', '^$[]')
CHR_IDX, POSITION_IDX, REF_NUC_IDX, N_READS_1_IDX, NUC_1_IDX, N_READS_2_IDX, NUC_2_IDX = range(7)
REF_IDX, TAXA1_IDX, TAXA2_IDX = 0, 1, 2
FLANK = 2  # 2 bases on each side for 5-mer

def parse_args():
    parser = argparse.ArgumentParser(description="Extract 5-mer mutations from pileup file.")
    parser.add_argument("reference", help="Reference species name")
    parser.add_argument("taxa1", help="Taxa 1 name")
    parser.add_argument("taxa2", help="Taxa 2 name")
    parser.add_argument("--pileup-dir", required=True, help="Directory with pileup.gz file")
    parser.add_argument("--output-dir", required=True, help="Directory for JSON output")
    parser.add_argument("--no-cache", action="store_true", help="Force recomputation")
    return parser.parse_args()

def parse_line(line):
    parts = line.strip().split('\t')
    return parts[:5] + parts[6:-1] if len(parts) >= 9 else None

def all_same(seq):
    return len(seq) > 0 and all(b == seq[0] for b in seq)

def get_nuc(field):
    cleaned = field.translate(REMOVE_CHARS)
    return cleaned[0].upper() if cleaned else 'N'

def quality_check(fields):
    return fields and '*' not in fields[NUC_1_IDX] and '*' not in fields[NUC_2_IDX] and \
           all_same(fields[NUC_1_IDX].translate(REMOVE_CHARS)) and \
           all_same(fields[NUC_2_IDX].translate(REMOVE_CHARS))

def extract_5mer(fields_list):
    sequences = [[], [], []]
    for fields in fields_list:
        ref_nuc = get_nuc(fields[REF_NUC_IDX])
        nuc1 = get_nuc(fields[NUC_1_IDX])
        nuc2 = get_nuc(fields[NUC_2_IDX])
        nuc1 = ref_nuc if nuc1 in {',', '.'} else nuc1
        nuc2 = ref_nuc if nuc2 in {',', '.'} else nuc2
        sequences[REF_IDX].append(ref_nuc)
        sequences[TAXA1_IDX].append(nuc1)
        sequences[TAXA2_IDX].append(nuc2)
    return sequences

def detect_mutation_5mer(five_mers):
    t1_mut, t2_mut = None, None
    flank_indices = [0, 1, 3, 4]
    center = 2

    # Require flanks to match across all 3
    for i in flank_indices:
        if not (five_mers[REF_IDX][i] == five_mers[TAXA1_IDX][i] == five_mers[TAXA2_IDX][i]):
            return None, None

    ref_base = five_mers[REF_IDX][center]
    t1_base = five_mers[TAXA1_IDX][center]
    t2_base = five_mers[TAXA2_IDX][center]

    context = ''.join(five_mers[REF_IDX])

    if t1_base == ref_base and t2_base != ref_base:
        t2_mut = f"{context[:2]}[{ref_base}>{t2_base}]{context[3:]}"
    elif t2_base == ref_base and t1_base != ref_base:
        t1_mut = f"{context[:2]}[{ref_base}>{t1_base}]{context[3:]}"
    return t1_mut, t2_mut

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pileup_path = os.path.join(args.pileup_dir, f"{args.reference}__{args.taxa1}__{args.taxa2}.pileup.gz")
    json1_path = os.path.join(args.output_dir, f"{args.taxa1}__{args.taxa2}__{args.reference}__5mers.json")
    json2_path = os.path.join(args.output_dir, f"{args.taxa2}__{args.taxa1}__{args.reference}__5mers.json")

    if all(os.path.exists(p) for p in [json1_path, json2_path]) and not args.no_cache:
        print("5-mer mutation files exist. Skipping.")
        return

    species_mut1 = defaultdict(int)
    species_mut2 = defaultdict(int)


    with gzip.open(pileup_path, 'rt') as f:
        window = [None] * (2 * FLANK + 1)
        qc = [False] * (2 * FLANK + 1)
        for i in range(2 * FLANK):
            window[i] = parse_line(f.readline())
            qc[i] = quality_check(window[i])

        for line in f:
            window = window[1:] + [parse_line(line)]
            qc = qc[1:] + [quality_check(window[-1])]

            if all(qc):
                five_mers = extract_5mer(window)
                mut1, mut2 = detect_mutation_5mer(five_mers)
                chrom = window[FLANK][CHR_IDX]
                pos = int(window[FLANK][POSITION_IDX])
                ref_base = five_mers[REF_IDX][FLANK]

                if mut1:
                    species_mut1[mut1] += 1
                if mut2:
                    species_mut2[mut2] += 1
    
    with open(json1_path, 'w') as f:
        json.dump(species_mut1, f, indent=2)
    with open(json2_path, 'w') as f:
        json.dump(species_mut2, f, indent=2)
    print(f"Written: {json1_path}, {json2_path}")

if __name__ == "__main__":
    main()
