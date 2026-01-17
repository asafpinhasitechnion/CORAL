import sys
import os
import subprocess
from pathlib import Path
import json
from datetime import datetime

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
    if len(sys.argv) < 7:
        print("Usage: run_pipeline.py <taxa1_name> <taxa1_accession> <taxa2_name> <taxa2_accession> <outgroup_name> <outgroup_accession> [optional args...]")
        sys.exit(1)

    t1_name, t1_acc, t2_name, t2_acc, out_name, out_acc = sys.argv[1:7]
    optional_args = sys.argv[7:]

    # === Argument groups ===
    download_args, index_args = [], []
    align_filter_args, pileup_args, mutation_args, triplet_args, interval_args, normalize_mutations_args = [], [], [], [], [], []
    global GENOMIC_PLOTS
    GENOMIC_PLOTS = False
    suffix = ""

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
        elif arg == '--genomic-position-plots':
            GENOMIC_PLOTS = True
            i += 1
        elif arg == '--suffix':
            suffix = optional_args[i + 1]
            i += 2
        elif arg == '--divergence-time':
            normalize_mutations_args.extend([arg, optional_args[i + 1]])
            i += 2
        else:
            print(f"Unknown argument: {arg}")
            sys.exit(1)

    return {
        "t1_name": t1_name,
        "t1_acc": t1_acc,
        "t2_name": t2_name,
        "t2_acc": t2_acc,
        "out_name": out_name,
        "out_acc": out_acc,
        "download_args": download_args,
        "index_args": index_args,
        "align_filter_args": align_filter_args,
        "pileup_args": pileup_args,
        "mutation_args": mutation_args,
        "triplet_args": triplet_args,
        "interval_args": interval_args,
        "suffix": suffix,
        "normalize_mutations_args": normalize_mutations_args
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

    run_id = f"{args['t1_name']}__{args['t2_name']}__{args['out_name']}"
    if "suffix" in args:
        run_id += '_' + args["suffix"]

    base_output_dir = Path("../Output") / run_id
    base_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Run ID: {run_id}")
    print(f"Base output directory: {base_output_dir}")
    write_metadata(args, base_output_dir)

    # === GENOME DOWNLOADS ===
    print("Downloading genomes...")
    run_cmd(["bash", "download_genome.sh", args["t1_name"], args["t1_acc"], str(base_output_dir)] + args["download_args"])
    run_cmd(["bash", "download_genome.sh", args["t2_name"], args["t2_acc"], str(base_output_dir)] + args["download_args"])
    run_cmd(["bash", "download_genome.sh", args["out_name"], args["out_acc"], str(base_output_dir)] + args["download_args"])

    # REFERENCE INDEXING
    print(f"Indexing outgroup genome: {args['out_name']}")
    run_cmd(["bash", "index_reference_genome.sh", args["out_name"], str(base_output_dir)] + args["index_args"])

    # ALIGNMENTS
    print(f"Aligning {args['t1_name']} to {args['out_name']}")
    run_cmd(["bash", "customizable_align_and_filter.sh", args["t1_name"], args["out_name"], str(base_output_dir)] + args["align_filter_args"])

    print(f"Aligning {args['t2_name']} to {args['out_name']}")
    run_cmd(["bash", "customizable_align_and_filter.sh", args["t2_name"], args["out_name"], str(base_output_dir)] + args["align_filter_args"])

    print(f"Alignment and filtering complete for {run_id}")

    # PILEUP
    print("Creating pileup...")
    run_cmd(["bash", "create_pileup.sh", args["out_name"], args["t1_name"], args["t2_name"], str(base_output_dir)] + args["pileup_args"])

    # MUTATIONS
    print("Extracting mutations...")
    run_cmd(["python3", "extract_mutations.py", args["out_name"], args["t1_name"], args["t2_name"],
             "--pileup-dir", str(base_output_dir),
             "--output-dir", str(base_output_dir / "Mutations")] + args["mutation_args"])

    # TRIPLETS
    print("Extracting triplet contexts...")
    run_cmd(["python3", "extract_triplets.py", args["out_name"], args["t1_name"], args["t2_name"],
             "--pileup-dir", str(base_output_dir),
             "--output-dir", str(base_output_dir / "Triplets")] + args["triplet_args"])

    # NORMALIZATION
    print("Normalizing mutation spectra...")
    run_cmd(["python3", "normalize_extracted_mutations.py", "--input-dir", str(base_output_dir)] + args["normalize_mutations_args"])

    # PLOTTING
    print("Plotting mutation and triplet spectra...")
    run_cmd(["python3", "plot_spectra.py", "--input-dir", str(base_output_dir / "Tables")])

    # INTERVAL EXTRACTION
    print("Extracting genomic intervals...")
    interval_dir = base_output_dir / "Intervals"
    interval_dir.mkdir(exist_ok=True)
    if GENOMIC_PLOTS:
        run_cmd([
            "python3", "bam_to_intervals.py",
            str(base_output_dir / "BAMs" / f"{args['t1_name']}_to_{args['out_name']}.bam"),
            str(interval_dir)
        ] + args["interval_args"])

        run_cmd([
            "python3", "bam_to_intervals.py",
            str(base_output_dir / "BAMs" / f"{args['t2_name']}_to_{args['out_name']}.bam"),
            str(interval_dir)
        ] + args["interval_args"])


        # === FAI FILE LOCATION ===
        fai_file = base_output_dir / args["out_name"] / f"{args['out_name']}.fasta.fai"
        if not fai_file.exists():
            print(f"FAI file not found: {fai_file}")
            sys.exit(1)

        top_chroms = get_top_n_chromosomes(fai_file, n=5)
        print("Plotting coverage and mutation density for top chromosomes...")
        for chrom in top_chroms:
            print(f"Plotting for {chrom}...")

            # Coverage plot
            run_cmd([
                "python3", "plot_genomic_coverage.py",
                "--interval-dir", str(interval_dir),
                "--fai-file", str(fai_file),
                "--chromosome", chrom,
                "--output", str(base_output_dir / f"Plots/coverage_{chrom}.png")
            ])

            # Mutation density plot: all mutations
            run_cmd([
                "python3", "plot_genomic_mutation_distribution.py",
                "--mutation-dir", str(base_output_dir / "Mutations"),
                "--coverage-dir", str(interval_dir),
                "--fai-file", str(fai_file),
                "--chromosome", chrom,
                "--output-dir", str(base_output_dir / "Plots")
            ])

            # Mutation density plot: C>T in [ACTG]CG context
            run_cmd([
                "python3", "plot_genomic_mutation_distribution.py",
                "--mutation-dir", str(base_output_dir / "Mutations"),
                "--coverage-dir", str(interval_dir),
                "--fai-file", str(fai_file),
                "--chromosome", chrom,
                "--mutation-category", r"[ACTG][C>T]G",
                "--output-dir", str(base_output_dir / "Plots")
            ])

                    # === CLEANUP ===
        if "--remove-temp" in sys.argv:
            print("Cleaning up intermediate files...")
            run_cmd([
                "bash", "cleanup_output.sh", str(base_output_dir),
                "--pileup", "--genomes",
                "--t1", args["t1_name"],
                "--t2", args["t2_name"],
                "--out", args["out_name"]
            ])


if __name__ == "__main__":
    main()
