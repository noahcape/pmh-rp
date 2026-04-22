#!/usr/bin/env python3

import subprocess
from pathlib import Path
import re

# ---------------- CONFIG ----------------

working_dir = "/n/fs/ragr-research/projects/pmh-rp/clone_trees"
clone_tree_tlp = (
    "/n/fs/ragr-research/projects/pmh-rp/scripts/dag_resolution_labeling.py"
)
metient = (
    "/n/fs/ragr-research/projects/pmh-rp/scripts/processing/run_metient_clone_trees.py"
)
process_metient = (
    "/n/fs/ragr-research/projects/pmh-rp/scripts/processing/process_metient_output.py"
)
draw_graph = "/n/fs/ragr-research/projects/pmh-rp/scripts/plots/draw_colored_tree.py"

data = ["A22", "A10", "CRUK0063", "P3"]
data = ["CRUK0063"]
# write the correct root per data then interleave when solving
roots = ["prostate", "prostate", "p.lung", "ROv"]
roots = ["p.lung"]


# helpers
def run_and_time(cmd, outfile):
    full_cmd = ["/usr/bin/time", "-v"] + cmd
    with open(outfile, "w") as f:
        subprocess.run(full_cmd, stdout=f, stderr=subprocess.STDOUT, check=True)


def mkdir_exist_ok(path):
    path.mkdir(parents=True, exist_ok=True)
    return path


def convert_metient_tree(txt_path, tsv_path):
    with open(txt_path) as fin, open(tsv_path, "w") as fout:
        for line in fin:
            parts = line.strip().split()
            if len(parts) >= 2:
                fout.write(f"s{parts[0]}\ts{parts[1]}\n")


def draw_migration_graph(tree, labeling, out, format: str, fout, tabs=False):
    with open(fout, "a") as f:
        proc = [
            "mamba",
            "run",
            "-n",
            "tlp",
            "python",
            str(draw_graph),
            str(tree),
            str(labeling),
            "-o",
            str(out),
            "-f",
            format,
            "-m",
        ]

        if tabs:
            proc = proc + ["-t"]

        subprocess.run(proc, stdout=f, stderr=subprocess.STDOUT, check=True)


# Running each method
def run_tlp(tree, labeling, root, tlp_out, timing_file):
    """Run TLP clone tree method."""
    print("-- Running tlp ---")
    run_and_time(
        [
            "mamba",
            "run",
            "-n",
            "tlp",
            "python",
            str(clone_tree_tlp),
            str(tree),
            str(labeling),
            "-o",
            str(tlp_out),
            "-r",
            str(root),
        ],
        timing_file,
    )

    labeling_file = f"{tlp_out}_vertex_labeling.csv"
    edge_file = f"{tlp_out}_apr.edgelist"
    draw_migration_graph(edge_file, labeling_file, tlp_out, "edgelist", timing_file)


def run_metient(tree, labeling, root, data_dir, metient_dir, timing_file):
    print("-- Running metient ---")
    """Run Metient pipeline (main + postprocess)."""
    # main run
    run_and_time(
        [
            "mamba",
            "run",
            "-n",
            "metient_gpu",
            "python",
            str(metient),
            str(tree),
            str(labeling),
            "-x",
            "1",
            "-r",
            str(root),
            "-m",
            str(metient_dir),
            "-o",
            str(metient_dir),
        ],
        timing_file,
    )

    # convert tree
    txt_tree = metient_dir / "metient_tree.txt"
    tsv_tree = metient_dir / "metient_tree.tsv"
    convert_metient_tree(txt_tree, tsv_tree)

    # postprocess
    with open(timing_file, "a") as f:
        subprocess.run(
            [
                "python",
                str(process_metient),
                str(metient_dir),
                str(timing_file),
                str(tsv_tree),
                str(tree),
                "empty",
            ],
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
        )

    # draw results
    pattern = re.compile(r"(\d+)_metient_edgelist\.tsv")

    for edge_file in metient_dir.glob("*_metient_edgelist.tsv"):
        match = pattern.match(edge_file.name)

        if not match:
            continue

        idx = match.group(1)
        labeling_file = f"{idx}_metient_labeling.csv"
        edge_file = Path(edge_file).name
        draw_migration_graph(
            metient_dir / edge_file,
            metient_dir / labeling_file,
            metient_dir / f"{idx}",
            "edgelist",
            timing_file,
        )


def run_mach2(tree, labeling, root, mach2_dir, timing_file):
    print("-- Running mach2 ---")
    """Run MACH2 solver."""
    run_and_time(
        [
            "mamba",
            "run",
            "-n",
            "mach2",
            "mach2",
            str(tree),
            str(labeling),
            "-o",
            str(mach2_dir),
            "-p",
            str(root),
            "--max_solutions",
            "10",
        ],
        timing_file,
    )

    # draw results
    pattern = re.compile(r"([a-zA-Z].+)-T-(\d+).refined\.tree")

    for edge_file in mach2_dir.glob("*.refined.tree"):
        match = pattern.match(edge_file.name)

        if not match:
            continue

        root_label = match.group(1)
        idx = match.group(2)
        labeling_file = f"{root_label}-T-{idx}.location.labeling"
        edge_file = Path(edge_file).name
        draw_migration_graph(
            mach2_dir / edge_file,
            mach2_dir / labeling_file,
            mach2_dir / f"{idx}",
            "edgelist",
            timing_file,
            tabs=True,
        )


def main():
    for d, root in zip(data, roots):
        data_dir = Path(working_dir) / d
        tree = data_dir / f"{d}.tree"
        labeling = data_dir / f"{d}.observed.labeling"
        mach2_dir = mkdir_exist_ok(data_dir / "mach2")
        metient_dir = mkdir_exist_ok(data_dir / "metient")
        tlp_dir = mkdir_exist_ok(data_dir / "tlp")
        tlp_out = tlp_dir / "tlp_output"

        mach2_timing_file = mach2_dir / "mach2_timing.txt"
        tlp_timing_file = tlp_dir / "tlp_clone_tree_timing.txt"
        metient_timing_file = metient_dir / "metient_timing.txt"

        # run_tlp(tree, labeling, root, tlp_out, tlp_timing_file)
        # run_metient(tree, labeling, root, data_dir, metient_dir, metient_timing_file)
        run_mach2(tree, labeling, root, mach2_dir, mach2_timing_file)

    print("Done.")


if __name__ == "__main__":
    main()
