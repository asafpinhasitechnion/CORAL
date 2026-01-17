import os
import json
import argparse
from collections import defaultdict
import pandas as pd
import re

# === Constants ===
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
valid_bases = {'A', 'C', 'G', 'T'}
mutation_pattern = re.compile(r"^[ACGT]\[[ACGT]>[ACGT]\][ACGT]$")

# === Functions ===
def load_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def save_json(obj, path):
    with open(path, 'w') as f:
        json.dump(obj, f, indent=2)

def collapse_triplets(triplet_dict):
    collapsed = defaultdict(int)
    for triplet, count in triplet_dict.items():
        if triplet[1] in {'G', 'A'}:
            rc = ''.join(complement.get(nuc, nuc) for nuc in reversed(triplet))
            collapsed[rc] += count
        else:
            collapsed[triplet] += count
    return collapsed

def get_complement(mutation):
    """Get the complement mutation."""
    comp_mutation = [complement[nuc] if nuc in complement else nuc for nuc in mutation]
    comp_mutation[0], comp_mutation[-1] = comp_mutation[-1], comp_mutation[0]
    return ''.join(comp_mutation)

def collapse_mutations(mutation_dict):
    """Collapse the mutations properly to reach 96 mutational categories."""
    collapsed_mutations = defaultdict(int)
    for mutation, count in mutation_dict.items():
        if mutation[2] in {'A', 'G'}:
            collapsed_mutations[get_complement(mutation)] += int(count)
        else:
            collapsed_mutations[mutation] += int(count)
    return collapsed_mutations

def filter_mutations_dict(d):
    return {k: v for k, v in d.items() if mutation_pattern.match(k)}

def filter_triplets_dict(d):
    return {k: v for k, v in d.items() if 'N' not in k and all(x in valid_bases for x in k)}

def normalize_by_triplets(mutations, triplets):
    return {
        k: v / triplets.get(f"{k[0]}{k[2]}{k[-1]}", 1) if triplets.get(f"{k[0]}{k[2]}{k[-1]}", 0) > 0 else 0
        for k, v in mutations.items()
    }

def scale_counts(d, target_sum=10000):
    total = sum(d.values())
    if total == 0:
        return {k: 0 for k in d}
    return {k: round(v / total * target_sum) for k, v in d.items()}

def main():
    parser = argparse.ArgumentParser(description="Normalize mutation spectra using triplet contexts.")
    parser.add_argument("--input-dir", required=True, help="Path to the run output directory (containing Mutations/ and Triplets/)")
    parser.add_argument("--output-dir", help="Optional output dir")
    parser.add_argument("--divergence-time", type=float, help="The divergence time of taxa1 and taxa2 (in million years ago, MYA)")

    args = parser.parse_args()

    mutation_dir = os.path.join(args.input_dir, "Mutations")
    triplet_dir = os.path.join(args.input_dir, "Triplets")
    output_dir = args.output_dir or os.path.join(args.input_dir, "Tables")

    # mutation_dir = '../Output/Plasmodium_fragile__Plasmodium_coatneyi__Plasmodium_gonderi/Mutations'
    # triplet_dir = '../Output/Plasmodium_fragile__Plasmodium_coatneyi__Plasmodium_gonderi/Triplets'
    # output_dir = '../Output/Plasmodium_fragile__Plasmodium_coatneyi__Plasmodium_gonderi/Tables'
    os.makedirs(output_dir, exist_ok=True)

    all_collapsed_mut = {}
    all_norm_mut = {}
    all_scaled_mut = {}
    all_triplets = {}

    for file in os.listdir(mutation_dir):
        if not file.endswith(".json"):
            continue

        mutation_path = os.path.join(mutation_dir, file)
        triplet_path = os.path.join(triplet_dir, file.replace("mutations", "triplets"))
        
        if not os.path.exists(triplet_path):
            continue

        mutations = filter_mutations_dict(load_json(mutation_path))
        triplets = filter_triplets_dict(load_json(triplet_path))

        collapsed_mut = collapse_mutations(mutations)
        collapsed_tri = collapse_triplets(triplets)

        all_collapsed_mut[file.split('.')[0]] = collapsed_mut
        all_triplets[file.split('.')[0]] = collapsed_tri

        norm_mut = normalize_by_triplets(collapsed_mut, collapsed_tri)
        scaled_mut = scale_counts(collapsed_mut)
        scaled_norm_mut = scale_counts(norm_mut)

        all_norm_mut[file.split('.')[0]] = scaled_norm_mut
        all_scaled_mut[file.split('.')[0]] = scaled_mut

    # Output combined 
    collapsed_mutations_df = pd.DataFrame(all_collapsed_mut)
    collapsed_mutations_df.to_csv(os.path.join(output_dir, "collapsed_mutations.tsv"), sep='\t')
    pd.DataFrame(all_norm_mut).to_csv(os.path.join(output_dir, "normalized_scaled.tsv"), sep='\t')
    pd.DataFrame(all_scaled_mut).to_csv(os.path.join(output_dir, "scaled_raw.tsv"), sep='\t')
    triplets_df = pd.DataFrame(all_triplets)
    triplets_df.to_csv(os.path.join(output_dir, "triplets.tsv"), sep='\t')
    
    mutations_per_triplet = collapsed_mutations_df.sum() / triplets_df.sum()
    print("Mutations per triplet from divergence:")
    for col in mutations_per_triplet.index:
        print(f"{col.split('__')[0]}: {mutations_per_triplet[col]:.2e}")

    if args.divergence_time:
        mutation_rates = mutations_per_triplet / (float(args.divergence_time) * 1_000_000)
        print("Estimated mutation rates per site per year:")
        for col in mutation_rates.index:
            print(f"{col.split('__')[0]}: {mutation_rates[col]:.2e}")



if __name__ == "__main__":
    main()
