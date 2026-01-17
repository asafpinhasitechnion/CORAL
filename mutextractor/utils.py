import subprocess
import sys

def log(message, verbose=True):
    if verbose:
        print(message, flush=True)


def run_cmd(cmd, shell=False, verbose = True):
    log(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}", verbose)
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        print(f"Command failed: {cmd}")
        sys.exit(result.returncode)

def get_top_n_chromosomes(fai_path, n=2):
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            chroms.append((fields[0], int(fields[1])))
    chroms.sort(key=lambda x: -x[1])
    return [c[0] for c in chroms[:n]]
