from ete3 import Tree
from Bio import Phylo
from io import StringIO
import sys
import os
import json

def parse_species_accession_from_newick(newick_str):
    tree = Phylo.read(StringIO(newick_str), "newick")
    species_accession_dict = {}
    for leaf in tree.get_terminals():
        # Expecting leaf names in the format: species|accession
        if "|" in leaf.name:
            species, accession = leaf.name.split("|", 1)
            species_accession_dict[species] = accession
        else:
            print(f"Leaf name '{leaf.name}' does not contain a '|' separator.")
            sys.exit(1)

    # Extract outgroup from the tree, or raise an error if no single outgroup exists
    root_clades = tree.root.clades
    outgroup = None
    for clade in root_clades:
        terminals = clade.get_terminals()
        if len(terminals) == 1:
            outgroup = terminals[0].name.split("|", 1)[0]
            break

    if outgroup is None:
        print("Could not determine a single outgroup from the Newick tree. Please ensure the tree is rooted and has a single outgroup.")
        sys.exit(1)

    return species_accession_dict, outgroup




def annotate_tree_with_indices(newick_str, outgroup_name, file_path=None):
    """
    Traverses the Newick tree using ete3, annotates each node with an index.
    Terminal nodes: 0 - outgroup, 1, 2, ... for the rest.
    Internal nodes: Node(idx).
    Returns:
      - tree (with .index and .custom_name on each node)
      - terminal_mapping: {species_name <-> idx} (bidirectional for terminals only)
    If file_path is provided, saves the annotated tree (as newick) and the mapping (as json)
    with suffixes '_annotated.nwk' and '_mapping.json'.
    """
    # Step 1: Load tree
    tree = Tree(newick_str, format=1)

    # Step 2: Clean names (strip suffix after "|")
    for leaf in tree.iter_leaves():
        if '|' in leaf.name:
            leaf.name = leaf.name.split('|', 1)[0]

    # Step 3: Sort terminals (outgroup first)
    terminals = tree.get_leaves()
    sorted_terminals = [t for t in terminals if t.name == outgroup_name] + \
                       sorted([t for t in terminals if t.name != outgroup_name], key=lambda x: x.name)

    # Step 4: Assign terminal indices and build mapping
    terminal_mapping = {}
    for idx, node in enumerate(sorted_terminals):
        node.add_feature("index", idx)
        node.add_feature("custom_name", node.name)
        species_name = node.name
        terminal_mapping[idx] = species_name
        terminal_mapping[species_name] = idx

    # Step 5: Assign internal node indices
    next_internal_idx = len(sorted_terminals)
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("index", next_internal_idx)
            node.add_feature("custom_name", f"Node({next_internal_idx})")
            next_internal_idx += 1

    # Step 6: Output to files (if needed)
    if file_path is not None:
        # Temporarily rename nodes to custom_name for writing
        original_names = {}
        for node in tree.traverse():
            original_names[node] = node.name
            node.name = getattr(node, "custom_name", node.name)

        annotated_tree_path = f"{os.path.splitext(file_path)[0]}_annotated.nwk"
        tree.write(format=1, outfile=annotated_tree_path)

        # Restore original names
        for node in tree.traverse():
            node.name = original_names[node]

        # Save mapping as JSON
        mapping_path = f"{os.path.splitext(file_path)[0]}_mapping.json"
        with open(mapping_path, "w") as f:
            json.dump(terminal_mapping, f, indent=2)

        print(f"Annotated tree saved to {annotated_tree_path}")
        print(f"Terminal mapping saved to {mapping_path}")

    return tree, terminal_mapping

def save_annotated_tree(tree, path):
    """
    Save an ETE3 tree including internal node names (custom_name).
    """
    # Temporarily assign names for all nodes
    original_names = {}
    for node in tree.traverse():
        original_names[node] = node.name
        node.name = getattr(node, "custom_name", node.name)

    tree.write(format=1, outfile=path)

    # Restore original names
    for node in tree.traverse():
        node.name = original_names[node]