#!/usr/bin/env python3
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
import numpy as np
import re

plt.rcParams.update({
    "font.size": 22,          # base font size
    "axes.titlesize": 22,
    "axes.labelsize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 15,
    "figure.titlesize": 20
})


# Helpers
def mean_err(x):
    """Return (mean, err). If only one value, err=None to suppress error bar"""
    if isinstance(x, (list, tuple, np.ndarray)):
        arr = np.asarray(x, dtype=float)
        if len(arr) == 1:
            return arr[0], None
        return arr.mean(), arr.std()   # or arr.std()/sqrt(len(arr)) for SEM
    else:
        return float(x), None


# === CONFIG ===
BASE_DIR = "/n/fs/ragr-research/projects/pmh-rp/laml"
DIR_PATTERN = r"m5k_lg\d+"
LEAF_LABELING_FILE = "laml_leaf_labeling.csv"
LEAF_THRESHOLD = 0
MACH2_SUBDIR = "mach2"

# JSON files
REGULARIZED_FILE = "tlp_regularized_exact_results.json"
COMPARE_FILES = [
    "mach2_results.json",
    "metient_results.json",
    "tlp_normal_results.json",
    # "tlp_dag_results.json",
    # "tlp_tree_results.json"
]

DATA = ["m5k_lg68", "m5k_lg28", "m5k_lg30", "m5k_lg17", "m5k_lg20", "m5k_lg21"]
SKIP = ["m5k_lg29", "m5k_lg46", "m5k_lg78", "m5k_lg94"]


def process_mach2_dir(mach2_path):
    """Compute average migrations and edges across all .graph files in a mach2 directory"""
    migrations_list = []
    edges_list = []

    for file in os.listdir(mach2_path):
        if file.endswith(".graph"):
            file_path = os.path.join(mach2_path, file)
            try:
                df = pd.read_csv(file_path, sep="\t", header=None)
                migrations = df[2].sum()
                edges = len(df)
                migrations_list.append(int(migrations))
                edges_list.append(edges)
            except Exception as e:
                print(f"Failed to parse {file_path}: {e}")

    if not migrations_list:
        return None

    # Compute averages
    return {
        "migrations": migrations_list[0],
        "migration_pattern_edges": edges_list[0],
        "num_files": len(migrations_list)
    }

def process_metient_dir(metient_path):
    """Compute average migrations and edges across all .graph files in a mach2 directory"""
    migrations_list = []
    edges_list = []
    files = []

    for file in os.listdir(metient_path):
        if file.endswith("info.json"):
            file_path = os.path.join(metient_path, file)
            try:
                with open(file_path, "r") as f:
                    data = json.load(f)
                    migrations_list.append(data["migrations"])
                    edges_list.append(data["num_edges"])
                    files.append(file_path)
            except Exception as e:
                print(f"Failed to parse {file_path}: {e}")

    if not migrations_list:
        return None

    min_mig = migrations_list[0]
    min_edge = edges_list[0]
    file = None

    for (m, e, fp) in zip(migrations_list, edges_list, files):
        if (e + m) < (min_mig + min_edge):
            min_edge = e
            min_mig = m
            file = fp

    # Compute averages
    return {
        "migrations": min_mig,
        "migration_pattern_edges": min_edge,
        "num_files": len(migrations_list)
    }

def generate_json_for_all_dirs(base_dir, processing_fx, method):
    for d in os.listdir(base_dir):
        dir_path = os.path.join(base_dir, d)
        if not os.path.isdir(dir_path) or not re.match(DIR_PATTERN, d):
            continue

        path = os.path.join(dir_path, method)
        if os.path.exists(path) and os.path.isdir(path):
            summary = processing_fx(path)
            if summary:
                # Save JSON inside the directory, mimicking TLP naming
                output_file = os.path.join(dir_path, f"{method}_results.json")
                with open(output_file, "w") as f:
                    json.dump(summary, f, indent=4)

def parse_json_file(file_path):
    with open(file_path, "r") as f:
        return json.load(f)

def collect_data(base_dir):
    """Collect migrations, edges, and compute deltas vs regularized"""
    records = []

    for d in os.listdir(base_dir):
        dir_path = os.path.join(base_dir, d)
        if not os.path.isdir(dir_path) or not re.match(DIR_PATTERN, d):
            continue

        # Leaf filtering
        leaf_file = os.path.join(dir_path, LEAF_LABELING_FILE)
        if not os.path.exists(leaf_file):
            continue
        try:
            leaf_df = pd.read_csv(leaf_file)
        except Exception as e:
            print(f"Failed to read {leaf_file}: {e}")
            continue
        if len(leaf_df) <= LEAF_THRESHOLD:
            continue

        if d in SKIP:
            continue

        # Parse regularized
        reg_path = os.path.join(dir_path, REGULARIZED_FILE)
        if not os.path.exists(reg_path):
            continue
        try:
            reg_data = parse_json_file(reg_path)
            reg_total = reg_data.get("migrations", 0) + reg_data.get("migration_pattern_edges", 0)
            # Add regularized as a record
            records.append({
                "directory": d,
                "method": "regularized",
                "leaves": len(leaf_df),
                "migrations": reg_data.get("migrations", 0),
                "migration_pattern_edges": reg_data.get("migration_pattern_edges", 0),
                "delta_vs_regularized": 0  # zero vs itself
            })
        except Exception as e:
            print(f"Failed to parse {reg_path}: {e}")
            continue

        # Compare to other JSON files
        for comp_file in COMPARE_FILES:
            comp_path = os.path.join(dir_path, comp_file)
            if os.path.exists(comp_path):
                try:
                    comp_data = parse_json_file(comp_path)
                    migrations = comp_data.get("migrations", 0)
                    edges = comp_data.get("migration_pattern_edges", 0)
                    comp_total = migrations + edges
                    delta = comp_total - reg_total
                    
                    records.append({
                        "directory": d,
                        "method": comp_file.replace("_results.json", ""),
                        "leaves": len(leaf_df),
                        "migrations": comp_data.get("migrations", 0),
                        "migration_pattern_edges": comp_data.get("migration_pattern_edges", 0),
                        "delta_vs_regularized": delta
                    })
                except Exception as e:
                    print(f"Failed to parse (here) {comp_path}: {e}")

    return pd.DataFrame(records)


def grouped_barplot(ax, df, ycol, directories, methods, palette, ylabel, title, truncate=None):
    width = 0.12
    x = np.arange(len(directories))

    for i, method in enumerate(methods):
        means = []
        errs = []

        for d in directories:
            row = df[(df["directory"] == d) & (df["method"] == method)]
            if len(row) == 0:
                means.append(np.nan)
                errs.append(None)
                continue

            val = row.iloc[0][ycol]
            m, e = mean_err(val)

            if truncate is not None:
                m = max(m, truncate)

            means.append(m)
            errs.append(e)

        # Convert errs so matplotlib ignores Nones
        yerr = [e if e is not None else 0 for e in errs]
        error_flags = [e is not None for e in errs]

        bars = ax.bar(
            x + i*width,
            means,
            width,
            color=palette[method],
            label=method,
            capsize=4
        )

        # # Add error bars manually only where needed
        # for j, bar in enumerate(bars):
        #     if error_flags[j]:
        #         ax.errorbar(
        #             bar.get_x() + bar.get_width()/2,
        #             means[j],
        #             yerr=yerr[j],
        #             fmt="none",
        #             ecolor="black",
        #             capsize=4
        #         )

    matches = [re.search(r"m5k_lg(\d+)", d) for d in directories]
    directories = [f"CP{int(m.group(1))}" for m in matches]

    ax.set_xticks(x + width*(len(methods)-1)/2)
    ax.set_xticklabels(directories, rotation=45)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()

    # horizontal lines
    ax.yaxis.grid(True, linestyle='--', color='gray', alpha=0.5)
    # optional: make them behind the bars
    ax.set_axisbelow(True)

def latex_table_edges_migrations(df, directories, methods, caption, label):
    # Map m5k_lgXX → CPXX
    cp_names = []
    for d in directories:
        m = re.search(r"m5k_lg(\d+)", d)
        cp_names.append(f"CP{int(m.group(1))}")

    table = {}
    for d in directories:
        table[d] = {}
        for mth in methods:
            row = df[(df["directory"] == d) & (df["method"] == mth)]
            if len(row) == 0:
                table[d][mth] = (np.nan, np.nan)
            else:
                mig = row.iloc[0]["migrations"]
                edg = row.iloc[0]["migration_pattern_edges"]
                mig, _ = mean_err(mig)
                edg, _ = mean_err(edg)
                table[d][mth] = (edg, mig)

    latex = []
    latex.append("\\begin{table*}[t]")
    latex.append("\\centering")
    latex.append("\\small")

    col_spec = "l" + "rr"*len(methods)
    latex.append(f"\\begin{{tabular}}{{{col_spec}}}")
    latex.append("\\toprule")

    # Method header
    header1 = [" "]
    for m in methods:
        header1.append(f"\\multicolumn{{2}}{{c}}{{{m}}}")
    latex.append(" & ".join(header1) + " \\\\")

    # Subheader
    header2 = ["Dataset"]
    for _ in methods:
        header2 += ["Edges", "Migr."]
    latex.append(" & ".join(header2) + " \\\\")
    latex.append("\\midrule")

    # Data rows
    for d, name in zip(directories, cp_names):
        row = [name]
        for m in methods:
            edg, mig = table[d][m]
            if np.isnan(edg):
                row += ["--", "--"]
            else:
                row += [f"{edg:.1f}", f"{mig:.1f}"]
        latex.append(" & ".join(row) + " \\\\")

    latex.append("\\bottomrule")
    latex.append("\\end{tabular}")
    latex.append(f"\\caption{{{caption}}}")
    latex.append(f"\\label{{{label}}}")
    latex.append("\\end{table*}")

    return "\n".join(latex)

def scatter_regularized_vs_edges_with_migration_diff(df, compare_to: str):
    output = f"figures/regularized_vs_{compare_to}_edges_migdiff.png"

    reg = df[df["method"] == "regularized"].set_index("directory")
    other = df[df["method"] == compare_to].set_index("directory")

    dirs = sorted(set(reg.index) & set(other.index))

    reg_edges, other_edges = [], []
    reg_migs, other_migs = [], []
    labels = []

    for d in dirs:
        red = reg.loc[d]["migration_pattern_edges"]
        rm  = reg.loc[d]["migrations"]
        me  = other.loc[d]["migration_pattern_edges"]
        mm  = other.loc[d]["migrations"]

        red, _ = mean_err(red)
        rm,  _ = mean_err(rm)
        me,  _ = mean_err(me)
        mm,  _ = mean_err(mm)

        reg_edges.append(red)
        other_edges.append(me)
        reg_migs.append(rm)
        other_migs.append(mm)

        m = re.search(r"m5k_lg(\d+)", d)
        labels.append(f"CP{int(m.group(1))}")

    reg_edges = np.array(reg_edges)
    other_edges = np.array(other_edges)
    reg_migs = np.array(reg_migs)
    other_migs = np.array(other_migs)

    # --- migration encoding ---
    mig_diff = np.abs(reg_migs - other_migs)

    colors = np.array([
        "tab:blue" if rm > om else
        "tab:purple" if rm < om else
        "tab:green"  # same number of migrations
        for rm, om in zip(reg_migs, other_migs)
    ])

    # robust size normalization
    size_min, size_max = 20, 200
    sizes = size_min + (size_max - size_min) * (
        mig_diff / mig_diff.max()
    )


    fig, ax = plt.subplots(figsize=(6, 6))
    # plt.rc('font', family='serif', serif=['Computer Modern Roman'])
    plt.rc('text', usetex=True)

    ax.scatter(
        reg_edges,
        other_edges,
        s=sizes,
        c=colors,
        alpha=0.8
    )

    # y = x reference
    ax.plot([0, 25], [0, 25], linestyle="--", linewidth=1)

    ax.set_xlim(0, 25)
    ax.set_ylim(0, 25)

    # --- add faint minor grid lines at every 1 ---
    ax.set_xticks(np.arange(0, 26, 1), minor=True)
    ax.set_yticks(np.arange(0, 26, 1), minor=True)

    ax.grid(True, linestyle="--", alpha=0.5)       # major grid (existing)
    ax.grid(which='minor', linestyle='--', alpha=0.15)  # faint minor grid

    ax.set_xlabel("Regularized edges")
    ax.set_ylabel(f"{compare_to.upper()} edges")

    if compare_to == "mach2":
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=f'reg $>$ {compare_to.upper()}', 
                markerfacecolor='tab:blue', markersize=8),
            Line2D([0], [0], marker='o', color='w', label=f'reg $=$ {compare_to.upper()}', 
                markerfacecolor='tab:green', markersize=8)
        ]
    else:
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=f'reg $>$ {compare_to.upper()}', 
                markerfacecolor='tab:blue', markersize=8),
            Line2D([0], [0], marker='o', color='w', label=f'reg $<$ {compare_to.upper()}', 
                markerfacecolor='tab:purple', markersize=8),
            Line2D([0], [0], marker='o', color='w', label=f'reg $=$ {compare_to.upper()}', 
                markerfacecolor='tab:green', markersize=8)
        ]

    ax.legend(handles=legend_elements, loc='lower right', frameon=False, title="Migrations")
    

    ax.set_axisbelow(True)

    # --- annotate largest migration difference ---
    imax = np.argmax(mig_diff)

    x = reg_edges[imax]
    y = other_edges[imax]
    rm = reg_migs[imax]
    om = other_migs[imax]

    ax.annotate(
        f"reg = {int(rm)}, {compare_to.upper()} = {int(om)}",
        xy=(x, y),
        xytext=(-6, 6),
        textcoords="offset points",
        fontsize=12,
        ha="right",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7"),
        arrowprops=dict(arrowstyle="->", linewidth=0.8)
    )

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print("Saved:", output)
    plt.close()


def scatter_regularized_vs_edges_with_migration_diff_ax(
    ax, df, compare_to: str
):
    reg = df[df["method"] == "regularized"].set_index("directory")
    other = df[df["method"] == compare_to].set_index("directory")

    dirs = sorted(set(reg.index) & set(other.index))

    reg_edges, other_edges = [], []
    reg_migs, other_migs = [], []

    for d in dirs:
        red, _ = mean_err(reg.loc[d]["migration_pattern_edges"])
        rm,  _ = mean_err(reg.loc[d]["migrations"])
        me,  _ = mean_err(other.loc[d]["migration_pattern_edges"])
        mm,  _ = mean_err(other.loc[d]["migrations"])

        reg_edges.append(red)
        other_edges.append(me)
        reg_migs.append(rm)
        other_migs.append(mm)

    reg_edges = np.array(reg_edges)
    other_edges = np.array(other_edges)
    reg_migs = np.array(reg_migs)
    other_migs = np.array(other_migs)

    mig_diff = np.abs(reg_migs - other_migs)

    colors = np.array([
        "tab:blue" if rm > om else
        "tab:purple" if rm < om else
        "tab:green"
        for rm, om in zip(reg_migs, other_migs)
    ])

    size_min, size_max = 20, 200
    sizes = size_min + (size_max - size_min) * (mig_diff / mig_diff.max())

    ax.scatter(reg_edges, other_edges, s=sizes, c=colors, alpha=0.8)

    ax.plot([0, 25], [0, 25], linestyle="--", linewidth=1)
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 25)

    ax.set_xticks(np.arange(0, 26, 1), minor=True)
    ax.set_yticks(np.arange(0, 26, 1), minor=True)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.grid(which="minor", linestyle="--", alpha=0.15)

    ax.set_xlabel("Regularized edges")
    if compare_to == "metient":
        ax.set_ylabel(f"{compare_to.capitalize()} edges")
    else:
        ax.set_ylabel(f"{compare_to.upper()} edges")

    ax.set_axisbelow(True)

    return (reg_migs - other_migs), size_min, size_max


def two_panel_scatter_plot(method1, method2, df):
    # fig = plt.figure(figsize=(14, 6))
    # gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.3], figure=fig)

    # axes = [fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])]
    # legend_ax = fig.add_subplot(gs[0,2])
    # legend_ax.axis('off')  # we only use this for the legend
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
    plt.rc('text', usetex=True)

    mig_diff_1, smin, smax = scatter_regularized_vs_edges_with_migration_diff_ax(
        axes[0], df, method1
    )
    mig_diff_2, _, _ = scatter_regularized_vs_edges_with_migration_diff_ax(
        axes[1], df, method2
    )

    # panel labels
    axes[0].set_title("a)", fontweight="bold", loc="left", pad=10)
    axes[1].set_title("b)", fontweight="bold", loc="left", pad=10)

    # Combine differences from both panels
    all_diff = np.concatenate([mig_diff_1, mig_diff_2])

    # signed differences: other - reg
    # assume mig_diff arrays were absolute; if you actually have signed diffs,
    # skip this. For absolute differences, we only show size, color can be neutral
    delta_m = all_diff   # replace with signed if you have it

    # define color mapping
    colors = np.array([
        "tab:blue" if d > 0 else
        "tab:green" if d == 0 else
        "tab:purple"  # if negative, for completeness
        for d in delta_m
    ])

    # define sizes (robust normalization across both panels)
    size_min, size_max = smin, smax
    sizes = size_min + (size_max - size_min) * (delta_m / delta_m.max())

    # unique differences present
    unique_deltas = np.unique(delta_m.astype(int))

    # build legend handles (combined size + color)
    legend_handles = []
    for d in unique_deltas:
        # color based on sign
        color = "tab:blue" if d > 0 else "tab:green" if d == 0 else "tab:purple"
        # size proportional to magnitude
        size = size_min + (size_max - size_min) * (abs(d) / delta_m.max() if delta_m.max() > 0 else 0)
        
        legend_handles.append(
            Line2D(
                [0], [0],
                marker='o',
                linestyle='None',
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=np.sqrt(size),
                label=rf"{d}",
                alpha=0.8,
            )
        )
    leg = fig.legend(
        handles=legend_handles,
        loc='lower center',
        ncol=len(legend_handles),   # number of columns in legend
        frameon=True,
        edgecolor='black',
        facecolor='white',
        title=rf"Migration difference (other − reg)",  # multi-line title
        # fontsize=10,
        # title_fontsize=11,
        handletextpad=0.5,
        labelspacing=0.5,
        handlelength=1.5,
        borderpad=0.5,
        borderaxespad=0.5
    )

    # left-align all labels
    for text in leg.get_texts():
        text.set_ha("left")
        text.set_fontweight('bold')

    # make title bold
    leg.get_title().set_fontweight('bold')


    plt.tight_layout(rect=[0, 0.18, 1, 1]) 
    plt.savefig("figures/regularized_vs_edges_combined.pdf", dpi=300)
    print("Saved to: figures/regularized_vs_edges_combined.pdf")
    plt.close()


def plot_all(df, filter=True):
    if df.empty:
        print("No data collected!")
        return
    
    if filter:    
        # Pivot only to decide which directories to keep
        wide = (
            df[df["method"].isin(["regularized", "tlp_normal"])]
            .pivot(index="directory", columns="method", values="migrations")
        )

        # Directories where regularized > tlp_normal
        keep_dirs = wide.index[wide["regularized"] > wide["tlp_normal"]]

        # Keep ALL methods for those directories
        df = df[df["directory"].isin(keep_dirs)]

    directories = df["directory"].unique()

    print(directories)
    methods = ["tlp_normal", "tlp_dag", "tlp_tree", "mach2", "metient", "regularized"]

    cmap = plt.get_cmap("tab10")
    palette = {m: cmap(i % 10) for i, m in enumerate(methods)}

    # palette = {
    #     "regularized": "#1f77b4",  # blue
    #     "tlp_dag": "#ff7f0e",       # orange
    #     "tlp_tree": "#2ca02c",      # green
    #     "tlp_normal": "#d62728",    # red
    #     "mach2": "#9467bd",         # purple
    #     "metient": "#17becf",       # cyan / teal
    # }

    fig, axes = plt.subplots(3, 1, figsize=(12, 15))

    # --- Migrations ---
    grouped_barplot(
        axes[0], df, "migrations",
        directories, methods, palette,
        ylabel="Migrations",
        title="Number of Migrations per Clone"
    )

    # --- Migration pattern edges ---
    grouped_barplot(
        axes[1], df, "migration_pattern_edges",
        directories, methods, palette,
        ylabel="Edges",
        title="Number of Migration Pattern Edges per Clone"
    )

    # --- Delta vs regularized ---
    delta_df = df[df["method"] != "regularized"]

    grouped_barplot(
        axes[2], delta_df, "delta_vs_regularized",
        directories,
        ["tlp_dag", "tlp_tree", "metient", "mach2", "tlp_normal"],
        palette,
        ylabel="Delta",
        title=f"Combined (migrations + edges) Difference to Regularized",
    )

    axes[2].axhline(0, color="black", linestyle="--", linewidth=1)

    plt.tight_layout()
    output_file = Path("lineage_tracing_data.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print("Saved:", output_file.resolve())
    plt.close()

if __name__ == "__main__":
    # generate_json_for_all_dirs(BASE_DIR, process_mach2_dir, "mach2")
    # generate_json_for_all_dirs(BASE_DIR, process_metient_dir, "metient")
    df = collect_data(BASE_DIR)

    directories = df["directory"].unique()
    # methods = ["tlp_normal", "tlp_dag", "tlp_tree", "mach2", "metient", "regularized"]
    methods = ["tlp_normal", "mach2", "metient", "regularized"]

    # print("\n=== EDGES + MIGRATIONS TABLE ===\n")
    # print(latex_table_edges_migrations(
    #     df,
    #     directories,
    #     methods,
    #     caption="Migration pattern edges and migrations per clone for each method.",
    #     label="tab:edges_migrations"
    # ))
    # print(df)
    # plot_all(df, False)

    # scatter_regularized_vs(df, "mach2")
    # scatter_regularized_vs(df, "metient")
    # scatter_regularized_vs_edges_with_migration_diff(df, "mach2")
    # scatter_regularized_vs_edges_with_migration_diff(df, "metient")

    two_panel_scatter_plot("mach2", "metient", df)
