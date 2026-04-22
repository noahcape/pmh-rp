#!/usr/bin/env python3
import subprocess
import json
from pathlib import Path
import csv
import os
from enum import auto, Enum

# ----------------------
# Configuration
# ----------------------
working_dir = Path("/n/fs/ragr-research/projects/pmh-rp")
sim_dir = Path("rust_sims_v2")  # directory of Rust simulations

tlp_script = working_dir / "scripts/tlp.py"
score_script = working_dir / "scripts/processing/score_result.py"
metient_script = working_dir / "scripts/processing/run_metient_modified.py"
process_metient_script = working_dir / "scripts/processing/process_metient_output.py"
out_base = working_dir / sim_dir

generations = [6,7,8,9,10]
# seeds = [1,2,3,4,5,6,7,8,9,10]
# generations = [6,7,8]
# seeds = [1,2,3,4,5]
seeds = [6,7,8,9,10]

class Method(Enum):
    TLP = auto(),
    MACH2 = auto(),
    METIENT = auto(),

# ----------------------
# Helpers
# ----------------------
def run_command(cmd, capture_output=False):
    print(f"Running: {' '.join(map(str, cmd))}")
    return subprocess.run(cmd, check=True, capture_output=capture_output, text=True)


def convert_edgelist(input_csv: Path, output_tsv: Path, round_digits: int = 0):
    """
    Converts a CSV edgelist (parent,child,length) to a TSV format with s<id> labels and rounded lengths.
    
    Example:
        0,1,1.423 -> 0 1 1
    """
    with input_csv.open("r", newline="") as f_in, output_tsv.open("w", newline="") as f_out:
        reader = csv.DictReader(f_in)
        for row in reader:
            parent = f"{row['parent']}"
            child = f"{row['child']}"
            length = round(float(row['length']), round_digits)
            f_out.write(f"{parent}\t{child}\t{length}\n")
            
def into_mach2_input_labeling(input_csv: Path, output_tsv: Path):
    """
    Convert CSV labeling file to MACH2 labeling format:
    - Remove header
    - Convert commas to tabs
    """
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with input_csv.open("r") as fin, output_tsv.open("w") as fout:
        next(fin)  # skip header
        for line in fin:
            fout.write(line.replace(",", "\t"))

def into_mach2_input_edgelist(input_tsv: Path, output_tsv: Path):
    """
    Extract first two columns of a TSV edgelist for MACH2.
    """
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with input_tsv.open("r") as fin, output_tsv.open("w") as fout:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            fout.write("\t".join(fields[:2]) + "\n")


def extract_stats(score_file: Path, stats_file: Path):
    if score_file.exists():
        with score_file.open() as f:
            data = json.load(f)
        stats = [
            data.get("inferred_migration_graph_num_edges"),
            data.get("inferred_parsimony_score"),
            data.get("pairwise_relations", {}).get("true_positives"),
            data.get("pairwise_relations", {}).get("true_negatives"),
            data.get("pairwise_relations", {}).get("false_positives"),
            data.get("pairwise_relations", {}).get("false_negatives"),
            data.get("pairwise_relations", {}).get("jaccard_index"),
            data.get("true_parsimony_score"),
            data.get("inferred_parsimony_score"),
            data.get("elapsed_time"),
        ]
        stats_file.parent.mkdir(parents=True, exist_ok=True)
        if not stats_file.exists():
            stats_file.write_text("num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips,time\n")
        with stats_file.open("a") as f:
            f.write(",".join(str(s) for s in stats) + "\n")

# ----------------------
# Helpers for method-specific directories and stats files
# ----------------------
def method_dir(sim_dir: Path, method: str) -> Path:
    """Directory for method output inside the simulation folder."""
    dir_path = sim_dir / method
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path

def stats_file_for_method(sim_dir: Path, method: str) -> Path:
    """CSV to collect results across seeds, placed in parent of sim_dir."""
    parent_dir = sim_dir.parent
    stats_path = parent_dir / f"{method}_results.csv"
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    if not stats_path.exists():
        stats_path.write_text("num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips,time\n")
    return stats_path

def tlp_type_f(regularize: bool, constraint: str) -> str:
    if regularize:
        return "tlp_reg"
    else:
        if constraint == "none":
            return "tlp_normal"
        elif constraint == "polyclonal_tree":
            return "tlp_tree"
        else:
            return "tlp_dag"

# ----------------------
# TLP Runner (with /usr/bin/time)
# ----------------------
def run_tlp(tree: Path, p_tree: Path, leaf_labeling: Path, full_labeling: Path,
            regularize: bool, constraint: str, sim_dir: Path):
    # method output directory
    tlp_type = tlp_type_f(regularize, constraint)
    out_dir = method_dir(sim_dir, tlp_type)
    out_prefix = out_dir / tlp_type
    
    timing_file = out_prefix.with_name(f"{out_prefix.name}_timing.txt")
    score_file = out_prefix.with_name(f"{out_prefix.name}_results.json")
    stats_file = stats_file_for_method(sim_dir, tlp_type)
    cmd = [
        "/usr/bin/time", "-v",
        "mamba", "run", "-n", "tlp", "python", str(tlp_script),
        "fast_machina",
        str(p_tree),
        str(leaf_labeling),
        "-o", str(out_prefix),
        "-e", str(int(regularize)),
        "-c", constraint
    ]
    with timing_file.open("w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)

    inferred_labeling = out_prefix.with_name(f"{out_prefix.name}_vertex_labeling.csv")

    cmd_score = [
        "mamba", "run", "-n", "tlp", "python", str(score_script),
        str(p_tree), str(tree),
        str(full_labeling),
        str(inferred_labeling),
        str(timing_file),
        "-o", str(score_file)
    ]
    run_command(cmd_score)
    extract_stats(score_file, stats_file)



# ----------------------
# MACH2 Runner (with /usr/bin/time)
# ----------------------
def run_mach2(tree: Path, p_tree: Path, leaf_labeling: Path, full_labeling: Path,
              sim_dir: Path, max_solutions: int = 10):
    """
    Runs MACH2 on a simulation, keeps outputs in sim_dir/mach2, captures timing, and
    aggregates all seeds into sim_dir.parent/mach2_results.csv
    """
    out_dir = method_dir(sim_dir, "mach2")
    out_prefix = out_dir / "mach2_run"
    
    timing_file = out_prefix.with_name(f"{out_prefix.name}_timing.txt")
    stats_file = stats_file_for_method(sim_dir, "mach2")
    
    # Convert CSV edgelist to MACH2 tree and labeling formats (stub)
    mach2_tree = out_prefix.with_name(f"{out_prefix.name}.tree")
    mach2_leaf = out_prefix.with_name(f"{out_prefix.name}.labeling")
    into_mach2_input_edgelist(p_tree, mach2_tree)
    into_mach2_input_labeling(leaf_labeling, mach2_leaf)
    
    cmd = [
        "/usr/bin/time", "-v",
        "mamba", "run", "-n", "mach2", "mach2",
        str(mach2_tree),
        str(mach2_leaf),
        "-o", str(out_dir),
        "--max_solutions", str(max_solutions)
    ]
    
    with timing_file.open("w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    
    # Score MACH2 outputs in TLP environment
    # MACH2 creates multiple solutions; loop through them
    for (f_index, file) in enumerate(os.listdir(out_dir)):
        if file.endswith(".location.labeling"):
            inferred_labeling = out_dir / file
            score_file = out_dir / f"{f_index}_scored.json"
            cmd_score = [
            "mamba", "run", "-n", "tlp", "python", str(score_script),
                str(mach2_tree),
                str(tree),
                str(full_labeling),
                str(inferred_labeling),
                str(timing_file),
                "-o", str(score_file)
            ]
            run_command(cmd_score)
            extract_stats(score_file, stats_file)


# ----------------------
# METIENT Runner (with /usr/bin/time)
# ----------------------
def run_metient(tree: Path, p_tree: Path, leaf_labeling: Path, full_labeling: Path,
                seed: int, sim_dir: Path):
    out_dir = method_dir(sim_dir, "metient")
    out_prefix = out_dir / f"metient_seed_{seed}"

    timing_file = out_prefix.with_name(f"{out_prefix.name}_timing.txt")
    stats_file = stats_file_for_method(sim_dir, "metient")
    metient_out = out_prefix
    metient_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "/usr/bin/time", "-v",
        "mamba", "run", "-n", "metient_gpu", "python", str(metient_script),
        str(p_tree),
        str(leaf_labeling),
        "-f", str(full_labeling),
        "-o", str(out_prefix),
        "-x", str(seed),
        "-r", "0",
        "-m", str(metient_out)
    ]
    with timing_file.open("w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)

    # process metient output
    cmd_process = [
        "mamba", "run", "-n", "metient_gpu", "python", str(process_metient_script),
        str(metient_out),
        str(timing_file),
        str(tree), # tree
        str(p_tree), # solved tree
        str(full_labeling),
        "-s", str(stats_file)
    ]
    run_command(cmd_process)
    
def run_sims(root_dir: Path, method: Method):
    for g in generations:
        for seed in seeds:
            path = root_dir / str(g) / str(seed)
            path_name = "sim"
            tree_in = path / f"{path_name}_edgelist.csv"
            tree_out = path / f"{path_name}_edgelist.tsv"
            convert_edgelist(tree_in, tree_out)
            p_tree_in = path / f"{path_name}_perturbed_edgelist.csv"
            p_tree_out = path / f"{path_name}_perturbed_edgelist.tsv"
            convert_edgelist(p_tree_in, p_tree_out)

            perturbed_leaf_labeling = path / f"{path_name}_perturbed_leaf_labeling.csv"
            full_labeling = path / f"{path_name}_vertex_labeling.csv"
            
            match method:
                case Method.TLP:
                    run_tlp(tree_out, p_tree_out, perturbed_leaf_labeling, full_labeling, regularize=False, constraint="none", sim_dir=path)
                    run_tlp(tree_out, p_tree_out, perturbed_leaf_labeling, full_labeling, regularize=True, constraint="none", sim_dir=path)
                    run_tlp(tree_out, p_tree_out, perturbed_leaf_labeling, full_labeling, regularize=False, constraint="polyclonal_tree", sim_dir=path)
                    run_tlp(tree_out, p_tree_out, perturbed_leaf_labeling, full_labeling, regularize=False, constraint="polyclonal_dag", sim_dir=path)
                case Method.MACH2:
                    run_mach2(tree_out, p_tree_out, perturbed_leaf_labeling, full_labeling, path)
                case Method.METIENT:
                    run_metient(tree_out, p_tree_out, perturbed_leaf_labeling, full_labeling, seed, path)
        
if __name__ == "__main__":
    sim_dir = working_dir / "rust_sims_v2"
    run_sims(sim_dir, Method.TLP)
    run_sims(sim_dir, Method.MACH2)
    run_sims(sim_dir, Method.METIENT)
    