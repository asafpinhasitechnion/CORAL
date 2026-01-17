import sys
import os
import subprocess
from pathlib import Path
import json
from datetime import datetime
from Bio import Phylo
from io import StringIO
from ete3 import Tree
from extract_multiple_species_mutations import extract_mutations
from multiple_species_utils import parse_species_accession_from_newick, annotate_tree_with_indices, save_annotated_tree
from run_phylip import run_phylip


def run_cmd(cmd, shell=False):
    print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        print(f"Command failed: {cmd}")
        sys.exit(result.returncode)


def write_metadata(args_dict, output_dir):
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "arguments": args_dict
    }

    metadata_path = Path(output_dir) / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"Metadata written to {metadata_path}")


def parse_args():
    # if len(sys.argv) < 3:
    #     print("Usage: run_pipeline.py <newick_species_tree> <run_id> [optional args...]")
    #     sys.exit(1)

    # newick_tree, run_id = sys.argv[1:3]
    # optional_args = sys.argv[3:]


    # newick_tree = "(Leptophobia_aripa|GCA_951799465.1, (Pieris_brassicae|GCF_905147105.1, (Pieris_napi|GCF_905475465.1, (Pieris_rapae|GCF_905147795.1, Pieris_mannii|GCA_028984075.1))));"
    # newick_tree = "(((Drosophila_miranda|GCF_003369915.1,Drosophila_pseudoobscura|GCF_009870125.1),Drosophila_helvetica|GCA_963969585.1),Drosophila_athabasca|GCA_008121215.1);"
    newick_tree = "(((Drosophila_sechellia|GCF_004382195.2,Drosophila_melanogaster|GCF_000001215.4),Drosophila_mauritiana|GCF_004382145.1),Drosophila_santomea|GCF_016746245.2);"

    run_id = 'drosophila2_run_mutiple_species'
    optional_args = ['--mapq', '1']

    # === Argument groups ===
    download_args, index_args = [], []
    align_filter_args, pileup_args, mutation_args, triplet_args, interval_args = [], [], [], [], []
    global GENOMIC_PLOTS
    GENOMIC_PLOTS = False
    reference = None

    i = 0
    while i < len(optional_args):
        arg = optional_args[i]
        if arg in ["--remove-temp"]:
            align_filter_args.append(arg)
            pileup_args.append(arg)
            i += 1
        elif arg == "--no-cache":
            # Add to all affected stages
            download_args.append(arg)
            index_args.append(arg)
            align_filter_args.append(arg)
            pileup_args.append(arg)
            mutation_args.append(arg)
            triplet_args.append(arg)
            interval_args.append(arg)
            i += 1
        elif arg == "--mapq":
            align_filter_args.extend([arg, optional_args[i + 1]])
            i += 2
        elif arg == "--no-full-mutations":
            mutation_args.append(arg)
            i += 1
        elif arg == "--aligner":
            align_filter_args.extend([arg, optional_args[i + 1]])
            i += 2
        elif arg == "--aligner-cmd":
            align_filter_args.extend([arg, optional_args[i + 1]])
            i += 2
        elif arg == '--reference':
            reference = optional_args[i + 1]
            i += 2
        elif arg == '--genomic-position-plots':
            GENOMIC_PLOTS = True
            i += 1
        else:
            print(f"Unknown argument: {arg}")
            sys.exit(1)

    return {
        "outgroup": reference,
        "newick_tree": newick_tree,
        "run_id": run_id,
        "download_args": download_args,
        "index_args": index_args,
        "align_filter_args": align_filter_args,
        "pileup_args": pileup_args,
        "mutation_args": mutation_args,
        "triplet_args": triplet_args,
        "interval_args": interval_args,
    }

def get_top_n_chromosomes(fai_path, n=2):
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            chroms.append((fields[0], int(fields[1])))
    chroms.sort(key=lambda x: -x[1])
    return [c[0] for c in chroms[:n]]



def main():
    args = parse_args()

    run_id = args['run_id']
    base_output_dir = Path("../Output") / run_id
    base_output_dir.mkdir(parents=True, exist_ok=True)
    newick_tree = args["newick_tree"]
    species_accession_dict, outgroup = parse_species_accession_from_newick(newick_tree)
    if args['outgroup']:
        outgroup = args['outgroup']

    print(f"Run ID: {run_id}")
    print(f"Base output directory: {base_output_dir}")
    write_metadata(args, base_output_dir)

    # === GENOME DOWNLOADS ===
    print("Downloading genomes...")
    for species, accession in species_accession_dict.items():
        run_cmd(["bash", "download_genome.sh", species, accession, str(base_output_dir)] + args["download_args"])

    # REFERENCE INDEXING
    print(f"Indexing outgroup genome: {outgroup}")
    run_cmd(["bash", "index_reference_genome.sh", outgroup, str(base_output_dir)] + args["index_args"])

    # ALIGNMENTS
    for species, accession in species_accession_dict.items():
        if species != outgroup:
            print(f"Aligning {species} to {outgroup}")
            run_cmd(["bash", "customizable_align_and_filter.sh", species, outgroup, str(base_output_dir)] + args["align_filter_args"])

    print(f"Alignment and filtering complete for {run_id}")

    tree, terminal_mapping = annotate_tree_with_indices(newick_tree, outgroup)
    # Extract terminal indices in order, skipping outgroup index 0
    ordered_taxa = [
        terminal_mapping[idx]
        for idx in sorted(k for k in terminal_mapping if isinstance(k, int) and k != 0)
    ]

    tree_path = os.path.join(base_output_dir, "annotated_tree.nwk")
    save_annotated_tree(tree, tree_path)
    with open(os.path.join(base_output_dir, "species_mapping.json"), 'w') as f:
        json.dump(terminal_mapping, f, indent=2)

    # PILEUP
    print("Creating pileup...")
    run_cmd(["bash", "create_multiple_species_pileup.sh", outgroup, str(base_output_dir), run_id] + ordered_taxa + args["pileup_args"])

    # # MUTATIONS
    print("Extracting mutations...")
    no_cache = "--no-cache" in args["mutation_args"]
    n_species = len(tree)
    pileup_file = os.path.join(base_output_dir, f"{run_id}.pileup.gz")
    extract_mutations(pileup_file, base_output_dir, n_species, tree, terminal_mapping, no_cache)
    # run_cmd(["python3", "extract_multiple_species_mutations.py", args["out_name"], args["t1_name"], args["t2_name"],
    #          "--pileup-dir", str(base_output_dir),
    #          "--output-dir", str(base_output_dir / "Mutations")] + args["mutation_args"])

    run_phylip(command='dnapars',
        df_path=  os.path.join(base_output_dir, "matching_bases.csv.gz"),
        tree_path=tree_path,
        output_dir=base_output_dir,
        prefix='multip_species_phylip',
        input_string='5\nY\n',
        mapping=terminal_mapping
    )
    
    

if __name__ == "__main__":
    main()
