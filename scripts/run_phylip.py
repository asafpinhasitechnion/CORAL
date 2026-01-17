import os
import shutil
import subprocess
import pandas as pd
from ete3 import Tree
from multiple_species_utils import annotate_tree_with_indices


def write_phylip_infile(df, outfile):
    irrelevant_cols = ['chromosome','position','left','right']
    taxa_cols = df.columns[~df.columns.isin(irrelevant_cols)]
    with open(outfile, 'w') as f:
        f.write(f"{len(taxa_cols)} {len(df)}\n")
        for col in taxa_cols:
            sequence = ''.join(df[col].astype(str).values).replace(' ', '-').upper()
            # f.write(f"{col.ljust(30)}{sequence}\n")
            f.write(f"{col}     {sequence}\n")


def write_intree(tree, outfile):
    original_names = {}
    for node in tree.traverse():
        original_names[node] = node.name
        node.name = getattr(node, "custom_name", node.name)
    tree.write(outfile=outfile, format=1)
    for node in tree.traverse():
        node.name = original_names[node]

import re

def extract_parsimony_score(outfile_path):
    """
    Extracts the parsimony score (total number of changes) from a PHYLIP outfile.

    Args:
        outfile_path (str): Path to the PHYLIP .outfile

    Returns:
        float: total number of changes
    """
    with open(outfile_path, 'r') as f:
        for line in f:
            if "requires a total of" in line:
                match = re.search(r"requires a total of\s+([0-9.]+)", line)
                if match:
                    return float(match.group(1))
    raise ValueError(f"Could not find parsimony score in {outfile_path}")

def run_phylip_command(df, output_dir, exe_path, tree=None, prefix="run1", phylip_input_args="Y\n", remove_infile=True):
    exe_path = os.path.abspath(exe_path)
    os.makedirs(output_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(output_dir)

    for fname in ["infile", "intree", "outfile", "outtree"]:
        if os.path.exists(fname):
            os.remove(fname)

    write_phylip_infile(df, "infile")
    print(f"PHYLIP input written to {os.path.abspath('infile')}")

    if tree:
        write_intree(tree, "intree")

    result = subprocess.run([exe_path], input=phylip_input_args, text=True, capture_output=True)

    with open("phylip_stdout.log", "w") as f:
        f.write(result.stdout)
    with open("phylip_stderr.log", "w") as f:
        f.write(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"PHYLIP run failed.\nExit code: {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    out_paths = {}
    if os.path.exists("outfile"):
        new_out = f"{prefix}.outfile"
        shutil.move("outfile", new_out)
        out_paths['outfile'] = os.path.abspath(new_out)
    if os.path.exists("outtree"):
        new_tree = f"{prefix}.outtree"
        shutil.move("outtree", new_tree)
        out_paths['outtree'] = os.path.abspath(new_tree)

    if remove_infile and os.path.exists("infile"):
        os.remove("infile")

    os.chdir(cwd)
    return out_paths


def run_phylip(command, df_path, tree_path, output_dir, prefix, input_string, mapping):

    df = pd.read_csv(df_path, index_col=0).astype(str)
    tree = Tree(tree_path, format=1) if tree_path else None

    if tree is not None and mapping is not None:
        for node in tree.iter_leaves():
            if node.name in mapping:
                node.name = f"taxa{mapping[node.name]}"
            else:
                raise ValueError(f"Species name '{node.name}' not found in mapping.")

    exe_path = os.path.abspath(f"./phylip-3.697/exe/{command}")

    no_tree_dir = os.path.join(output_dir, f"{prefix}_no_tree")
    tree_dir = os.path.join(output_dir, f"{prefix}_with_tree")

    print(f"Running PHYLIP {command} without a starting tree...")
    default_output = run_phylip_command(
        df,
        output_dir=no_tree_dir,
        exe_path=exe_path,
        tree=None,
        prefix=prefix,
        phylip_input_args=input_string
    )
    print("PHYLIP outputs (no tree):", default_output)

    print(f"Running PHYLIP {command} with a starting tree...")
    tree_output = run_phylip_command(
        df,
        output_dir=tree_dir,
        exe_path=exe_path,
        tree=tree,
        prefix="given_tree_run",
        phylip_input_args='U\n' + input_string
    )
    print("PHYLIP outputs (with tree):", tree_output)

    
    # === Compare scores ===
    default_score = extract_parsimony_score(default_output["outfile"])
    given_score = extract_parsimony_score(tree_output["outfile"])

    print("\nParsimony Score Comparison:")
    print(f" - Most Parsimonious Tree Score: {default_score}")
    print(f" - Given Tree Score: {given_score}")

    if given_score == default_score:
        print("The input tree is the most parsimonious.")
    elif given_score > default_score:
        ratio = given_score/default_score
        print(f"The input tree requires {ratio:.2f} times more changes than the most parsimonious tree.")
    else:
        print(f"Unexpected: input tree is more parsimonious than the optimal tree (check logic).")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run PHYLIP on mutation matrix.")
    parser.add_argument("--command", type=str, required=True,
                        help="PHYLIP command to run: 'dnapars', 'dnapenny', or 'dnamlk'")
    parser.add_argument("--input", type=str, default="../Output/test_run_mutiple_species/matching_bases.csv.gz",
                        help="Path to input CSV mutation matrix")
    parser.add_argument("--tree-input", type=str, default=None,
                        help="Path to input newick tree file")
    parser.add_argument("--output-dir", type=str, default="../Output/test_run_mutiple_species/phylip_run",
                        help="Directory to write PHYLIP outputs")
    parser.add_argument("--prefix", type=str, default="species_run",
                        help="Prefix for output files")
    parser.add_argument("--interactive", type=str, default="Y\n",
                        help="Interactive input for PHYLIP (default: just accept settings)")
    parser.add_argument("--outgroup", type=str, default=None,
                        help="Outgroup species name for tree annotation")
    args = parser.parse_args()
    _, mapping = annotate_tree_with_indices(args.tree_input, args.outgroup)

    run_phylip(
        command=args.command,
        df_path=args.input,
        tree_path=args.tree_input,
        output_dir=args.output_dir,
        prefix=args.prefix,
        input_string=args.interactive,
        mapping=mapping
        
    )
