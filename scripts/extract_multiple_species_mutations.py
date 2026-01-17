#!/usr/bin/env python3
from collections import defaultdict
import json
import os
import sys
import argparse
import gzip
import csv
import time
from ete3 import Tree
import pandas as pd
from multiple_species_utils import annotate_tree_with_indices, save_annotated_tree
from normalize_extracted_mutations import collapse_mutations
from plot_spectra import plot_mutations

# === Constants ===
REMOVE_CHARS = str.maketrans('', '', '^$[]')
CHR_IDX, POS_IDX, REF_NUC_IDX = 0, 1, 2
PREV_IDX, CUR_IDX, NEXT_IDX = 0, 1, 2

# === Argument Parser ===
def parse_args():
    parser = argparse.ArgumentParser(description="Extract positions with conserved flanks and differing middle bases across species.")
    parser.add_argument("--n-species", type=int, required=True, help="Number of species in the pileup file")
    parser.add_argument("--pileup-file", required=True, help="Path to .pileup.gz file")
    parser.add_argument("--output-dir", required=True, help="Path to output files")
    parser.add_argument("--newick-tree", required=True, help="Path to newick tree file")
    parser.add_argument("--outgroup", required=True, help="Name of outgroup species")
    parser.add_argument("--no-cache", action="store_true", help="Do not use cached mutations")

    # Optional hardcoded test mode (for debugging)
    if len(sys.argv) == 1:
        test_args = [
            
            "--pileup-file", "../Output/test_run_mutiple_species/Leptophobia_aripa__Pieris_brassicae__Pieris_mannii__Pieris_napi__Pieris_rapae.pileup.gz",
            "--n-species", "5",
            "--output-dir", "../Output/test_run_mutiple_species",
            "--newick-tree", "(Leptophobia_aripa, (Pieris_brassicae, (Pieris_napi, (Pieris_rapae, Pieris_mannii))));",
            "--outgroup", "Leptophobia_aripa"
        ]
        return parser.parse_args(test_args)
    return parser.parse_args()

# === Utilities ===
def all_same(seq):
    return len(seq) > 0 and all(ch == seq[0] for ch in seq)

def quality_check(fields):
    sample_fields = fields[3:]
    return (
        sample_fields
        and all('*' not in field for field in sample_fields)
        and all(all_same(field.translate(REMOVE_CHARS)) for field in sample_fields)
    )

def parse_line(line, n_species):
    parts = line.strip().split('\t')
    if len(parts) < n_species * 3:
        return None

    # Extract chromosome, position, reference base
    chrom, pos, ref_base = parts[:3]

    # Extract base strings: column 4 + 3*i for each sample
    base_calls = parts[4::3]

    # Normalize each base: take the first character or use ref if '.' or ','
    normalized = [
        base[0] if base and base[0] not in {',', '.'} else ref_base
        for base in base_calls
    ]
    return [chrom, pos, ref_base] + normalized

def detect_mutations(triplet_buffer):
    # Extract triplets from 3-line buffer (prev, curr, next)
    triplets = [fields[3:] for fields in triplet_buffer]

    prev_bases = triplets[PREV_IDX]
    curr_bases = triplets[CUR_IDX]
    next_bases = triplets[NEXT_IDX]

    if all_same(prev_bases) and all_same(next_bases) and len(set(curr_bases)) > 1:
        chrom = triplet_buffer[CUR_IDX][CHR_IDX]
        pos = triplet_buffer[CUR_IDX][POS_IDX]
        ref_base = triplet_buffer[CUR_IDX][REF_NUC_IDX]
        return [chrom, pos, prev_bases[0].upper(), next_bases[0].upper(), ref_base.upper()] + [base.upper() for base in curr_bases]
    return None



def recursive_state_check(node, row, mapping): 
    if node.is_leaf():
        node.add_feature("state", {row[f"taxa{mapping[node.name]}"]})
        return node.state

    # Recursively compute children's states
    left_state = recursive_state_check(node.children[0], row, mapping)
    right_state = recursive_state_check(node.children[1], row, mapping)

    intersect = left_state & right_state
    node_state = intersect if intersect else left_state | right_state
    node.add_feature("state", node_state)
    return node_state


def recursive_fitch(node, parent_state, row, mutation_dict, ambiguous_count, verbose=False):
    next_state = parent_state

    # Only if the parent state isn't compatible with the current state
    if parent_state not in node.state:
        if len(node.state) > 1:
            ambiguous_count += 1
            if verbose:
                print(f"{node.name} state is ambiguous: {node.state}. Cannot resolve mutations in this subtree.")
            return mutation_dict, 1
        else:
            next_state = list(node.state)[0]
            parent_name = node.up.custom_name if node.up else "ROOT"
            child_name = node.custom_name
            branch_key = f"{parent_name}→{child_name}"
            
            # Build mutation string
            mutation = f"{row['left']}[{parent_state}>{next_state}]{row['right']}"
            pos = row['position']
            chrom = row['chromosome']
            mutation_dict.setdefault(branch_key, []).append((chrom, pos, mutation))

    # Recurse on children if not a leaf
    if not node.is_leaf():
        mutation_dict, ambiguous_count = recursive_fitch(node.children[0], next_state, row, mutation_dict, ambiguous_count, verbose)
        mutation_dict, ambiguous_count = recursive_fitch(node.children[1], next_state, row, mutation_dict, ambiguous_count, verbose)

    return mutation_dict, ambiguous_count


def fitch(tree_root, row, mapping, mutation_dict, verbose=False):
    root_state = recursive_state_check(tree_root, row, mapping)
    ambiguous_count = 0
    if len(root_state) == 1:
        root_base = list(root_state)[0]
        mutation_dict, ambiguous_count = recursive_fitch(tree_root, root_base, row, mutation_dict, ambiguous_count, verbose)
    else:
        ambiguous_count += 1
        if verbose:
            print(f"Root state is ambiguous: {root_state}. Cannot resolve mutations.")

    return mutation_dict, ambiguous_count


def extract_mutations(pileup_file, output_dir, n_species, tree, mapping, no_cache=False):
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    plots_dir = os.path.join(output_dir, "Plots")
    os.makedirs(plots_dir, exist_ok=True)
    csv_dir = os.path.join(output_dir, "CSVs")
    os.makedirs(csv_dir, exist_ok=True)
    
    csv_output_path = os.path.join(output_dir, "matching_bases.csv.gz")
    data = []
    header = ["chromosome", "position", "left", "right"] + [f"taxa{i}" for i in range(n_species)]
    if os.path.exists(csv_output_path) and not no_cache:
        print(f"Loading cached mutations from {csv_output_path}")
        df = pd.read_csv(csv_output_path)
    else:
        with gzip.open(pileup_file, 'rt') as infile:
            # Initialize 3-line buffer for sliding window
            buffer = [None, parse_line(infile.readline(), n_species), parse_line(infile.readline(), n_species)]
            qc_flags = [False, quality_check(buffer[1]), quality_check(buffer[2])]

            for line in infile:
                buffer = [buffer[1], buffer[2], parse_line(line, n_species)]
                qc_flags = [qc_flags[1], qc_flags[2], quality_check(buffer[2])]

                if all(qc_flags):
                    result = detect_mutations(buffer)
                    if result:
                        data.append(result)

        df = pd.DataFrame(data, columns=header)
        df.to_csv(csv_output_path, index=False)
        print(f"Matching triplet mutation positions saved into {csv_output_path}")


    mutation_dict = defaultdict(list)
    ambiguous_counter = 0

    for _, row in df.iterrows():
        tree_copy = tree.copy()

        mutation_dict, ambiguous = fitch(
            tree_copy,
            row,
            mapping,
            mutation_dict,
            verbose=False
        )
        ambiguous_counter += ambiguous

    spectras_dicts = {}
    # === Step 3: Save mutations by branch ===
    for branch_key, mutations in mutation_dict.items():
        csv_path = os.path.join(csv_dir, f"{branch_key}.csv.gz")
        print(f"Saving mutations for {branch_key} to {csv_path}")
        df = pd.DataFrame(mutations, columns=["chromosome", "position", "mutation"])
        df.to_csv(csv_path, index=False, header=False, sep="\t", compression="gzip")
        mutation_spectra = dict(df['mutation'].value_counts())
        mutation_spectra = collapse_mutations(mutation_spectra)
        spectras_dicts[branch_key] = mutation_spectra
        spectra_plot_path = os.path.join(plots_dir, f"{branch_key}_spectra.png")
        print(f"Plotting mutation spectra for {branch_key} to {spectra_plot_path}")
        plot_mutations(pd.Series(mutation_spectra), spectra_plot_path, f"Mutation Spectra: {branch_key}")
    spectra_df = pd.DataFrame(spectras_dicts)
    print(f"Saving combined mutation spectras to {os.path.join(output_dir, 'mutation_spectras.tsv')}")
    spectra_df.to_csv(os.path.join(output_dir, "mutation_spectras.tsv"), sep="\t")


def memory_efficient_extract_mutations(pileup_file, output_dir, n_species, tree, mapping, no_cache=False):
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    csv_output_path = os.path.join(output_dir, "matching_bases.csv.gz")
    header = ["chromosome", "position", "left", "right"] + [f"taxa{i}" for i in range(n_species)]
    if not os.path.exists(csv_output_path) or no_cache:
        with gzip.open(csv_output_path, 'wt', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(header)
            with gzip.open(pileup_file, 'rt') as infile:
                # Initialize 3-line buffer for sliding window
                buffer = [None, parse_line(infile.readline(), n_species), parse_line(infile.readline(), n_species)]
                qc_flags = [False, quality_check(buffer[1]), quality_check(buffer[2])]

                for line in infile:
                    buffer = [buffer[1], buffer[2], parse_line(line, n_species)]
                    qc_flags = [qc_flags[1], qc_flags[2], quality_check(buffer[2])]

                    if all(qc_flags):
                        result = detect_mutations(buffer)
                        if result:
                            writer.writerow(result)
        print(f"Mutations written to {csv_output_path}")
    
    mutation_dict = {}
    ambiguous_counter = 0

    with gzip.open(csv_output_path, 'rt') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            mutation_update, ambiguous_count = fitch(tree.copy(), row, mapping, mutation_dict, verbose=False)
            ambiguous_counter += ambiguous_count

    print(f"Total ambiguous mutations: {ambiguous_counter}")

    # === Step 3: Save mutations by branch ===
    for branch_key, mutations in mutation_dict.items():
        output_path = os.path.join(output_dir, f"{branch_key}.csv.gz")
        print(f"Saving mutations for {branch_key} to {output_path}")
        df = pd.DataFrame(mutations, columns=["chromosome", "position", "mutation"])
        df.to_csv(output_path, index=False, header=False, sep="\t", compression="gzip")


# === Main Function ===
def main():
    args = parse_args()
    n_species = args.n_species
    pileup_file = args.pileup_file
    output_dir = args.output_dir
    newick_tree = args.newick_tree
    outgroup = args.outgroup
    no_cache = args.no_cache
    tree, mapping = annotate_tree_with_indices(newick_tree, outgroup)
    save_annotated_tree(tree, os.path.join(output_dir, "annotated_tree.nwk"))
    with open(os.path.join(output_dir, "species_mapping.json"), 'w') as f:
        json.dump(mapping, f, indent=2)
    extract_mutations(pileup_file, output_dir, n_species, tree, mapping, no_cache=no_cache)


if __name__ == "__main__":
    main()