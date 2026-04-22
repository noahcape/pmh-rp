import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
import numpy as np

matplotlib.use("Agg")

def parsimony_score_gap(dfs, labels, simulation_type, colors=None):
    # Compute percentage gap for each file and drop NaNs
    gaps = []
    for df in dfs:
        gap = df["ips"] - df["tps"]
        gaps.append(gap.dropna())
    
    # Compute mean gap for sorting by closeness to zero
    means = [g.mean() for g in gaps]
    medians = [g.median() for g in gaps]
    sorted_indices = np.argsort(np.abs(means))[::-1]  # furthest from zero first
    gaps = [gaps[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]
    
    means = [means[i] for i in sorted_indices]  # reorder means to match labels

    all_gaps = pd.concat(
        gaps,
        keys=labels,
        names=["Method", "Index"]
    )
    all_gaps = all_gaps.reset_index(level="Method")  # bring method name as column
    all_gaps.rename(columns={0: "Gap"}, inplace=True)
    all_gaps.to_csv("parsimony_gaps.csv", index=False)
    print("Saved all gaps to parsimony_gaps.csv")

    # Set colors
    n = len(gaps)
    if colors is None:
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i % 10) for i in range(n)]
    else:
        colors = [colors[i] for i in sorted_indices]
    # Ensure all colors are valid
    colors = [c if isinstance(c, (tuple, str)) else plt.get_cmap("tab10")(i % 10)
              for i, c in enumerate(colors)]

    # Create figure
    plt.figure(figsize=(8, 5))
    box = plt.boxplot(gaps, patch_artist=True, medianprops=dict(color='black'))

    # Color boxes
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Color whiskers and caps safely
    for i in range(n):
        # Whiskers
        box['whiskers'][2*i].set_color(colors[i])
        box['whiskers'][2*i+1].set_color(colors[i])
        box['whiskers'][2*i].set_linewidth(1.5)
        box['whiskers'][2*i+1].set_linewidth(1.5)
        # Caps
        box['caps'][2*i].set_color(colors[i])
        box['caps'][2*i+1].set_color(colors[i])
        box['caps'][2*i].set_linewidth(1.5)
        box['caps'][2*i+1].set_linewidth(1.5)
    # Medians
    for median in box['medians']:
        median.set_color('k')
        median.set_linewidth(2)

    plt.ylabel("Parsimony Score Gap (inferred - true)")
    plt.title(f"Parsimony Score Gap ({simulation_type})")
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    # Hide x-axis ticks
    plt.xticks([])

    # Create legend
    for label, color in zip(labels, colors):
        plt.plot([], [], color=color, linewidth=8, label=label)
    plt.legend(title="Method", loc="upper left")

    plt.tight_layout()
    plt.savefig("parsimony_score_gap.png", dpi=300, bbox_inches="tight")
    print("Saved:", Path("parsimony_score_gap.png").resolve())
    plt.close()


def f1_boxplot(dfs, labels, simulation_type, colors=None):
    def f1_score(df):
        tp = df["tp"]
        fp = df["fp"]
        fn = df["fn"]

        denom = 2 * tp + fp + fn
        return (2 * tp) / denom.where(denom != 0)

    # Compute F1 scores and drop NaNs
    f1_scores = [f1_score(df).dropna() for df in dfs]

    # Compute mean for sorting
    means = [s.mean() for s in f1_scores]
    sorted_indices = np.argsort(means)  # descending
    f1_scores = [f1_scores[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]
    
    # set colors
    n = len(f1_scores)
    if colors is None:
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i % 10) for i in range(n)]
    colors = [colors[i] for i in sorted_indices]  # reorder after sorting


    plt.figure(figsize=(8, 5))

    # Create boxplot
    box = plt.boxplot(f1_scores, patch_artist=True, boxprops=dict(facecolor='black'), medianprops=dict(color='black', linewidth=1))

    # Color each box
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    plt.ylabel("Migration Graph F1 Score")
    plt.ylim(-0.05, 1.05)
    plt.title(f"F1 Score Distribution Comparison ({simulation_type})")
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    # Hide x-axis labels
    plt.xticks([])

    # Create legend
    for label, color in zip(labels, colors):
        plt.plot([], [], color=color, label=label, linewidth=8)
    plt.legend(title="Method", loc="lower left")

    plt.tight_layout()
    plt.savefig("f1_score_boxplot.png", dpi=300, bbox_inches="tight")
    print("Saved:", Path("f1_score_boxplot.png").resolve())
    plt.close()

def f1_violinplot(dfs, labels, simulation_type, colors=None):
    def f1_score(df):
        tp = df["tp"]
        fp = df["fp"]
        fn = df["fn"]

        denom = 2 * tp + fp + fn
        return (2 * tp) / denom.where(denom != 0)

    # compute F1 scores and drop NaNs
    f1_scores = [f1_score(df).dropna() for df in dfs]

    # compute mean for each method and sort
    means = [s.mean() for s in f1_scores]
    sorted_indices = np.argsort(means)  # descending

    f1_scores = [f1_scores[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]

    n = len(f1_scores)
    x = np.arange(1, n + 1)

    if colors is None:
        # Use default colormap
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i % 10) for i in range(n)]
    colors = [colors[i] for i in sorted_indices]

    plt.figure(figsize=(8, 5))

    # Plot each violin individually with a different color
    for xi, scores, color in zip(x, f1_scores, colors):
        parts = plt.violinplot([scores], positions=[xi], showmeans=True, showmedians=True)
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        for pc in parts:
            if pc != "bodies":
                parts[pc].set_color('k')

    plt.ylabel("Migration Graph F1 Score")
    plt.ylim(-0.05, 1.05)
    plt.title(f"F1 Score Distribution Comparison ({simulation_type})")
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    # Hide x-axis ticks
    plt.xticks([])

    # Create legend
    for label, color in zip(labels, colors):
        plt.plot([], [], color=color, label=label, linewidth=8)
    plt.legend(title="Method", loc="lower left")

    plt.tight_layout()
    plt.savefig("f1_score_violin.png", dpi=300, bbox_inches="tight")
    print("Saved:", Path("f1_score_violin.png").resolve())
    plt.close()


migration_rate = [0.001, 0.0015, 0.002]
leaves = [50, 100, 500, 1000]
constraints = ["none", "polyclonal_dag", "polyclonal_tree"]

base = Path("/n/fs/ragr-research/projects/pmh-rp/simulations")
dfs = defaultdict(list)

def load_csv(path):
    return pd.read_csv(path).dropna(how="all")

for mrate in migration_rate:
    for num_leaves in leaves:
        for constraint in constraints:
            params = f"leaves={num_leaves}_constraint={constraint}_mrate={mrate}"
            root = base / params

            dfs["mach2"].append(load_csv(root / f"mach2_{constraint}.csv"))
            if num_leaves < 100:
                dfs["metient"].append(load_csv(root / "metient.csv"))
            dfs["tlp_regularized"].append(
                load_csv(root / f"tlp_{constraint}_r=1.csv")
            )

            # normal TLP
            if constraint == "none":
                dfs["tlp_normal"].append(
                    load_csv(root / "tlp_none_r=0.csv")
                )
            else:
                dfs["tlp_normal"].append(
                    load_csv(root / f"tlp_{constraint}.csv")
                )

            # constrainted TLP
            if constraint == "polyclonal_tree":
                dfs["tlp_tree"].append(
                    load_csv(root / f"tlp_{constraint}_r=0.csv")
                )
            elif constraint == "polyclonal_dag":
                dfs["tlp_dag"].append(
                    load_csv(root / f"tlp_{constraint}_r=0.csv")
                )

# Final concatenation
mach2 = pd.concat(dfs["mach2"], ignore_index=True)
metient = pd.concat(dfs["metient"], ignore_index=True)
tlp_regularized = pd.concat(dfs["tlp_regularized"], ignore_index=True)
tlp_normal = pd.concat(dfs["tlp_normal"], ignore_index=True)
tlp_tree = pd.concat(dfs["tlp_tree"], ignore_index=True)
tlp_dag = pd.concat(dfs["tlp_dag"], ignore_index=True)
constrained = pd.concat([tlp_tree, tlp_dag], ignore_index=True)

fname_labels = ["mach2", "metient", "regularized parsimony", "normal parsimony"]
fname_labels = ["mach2", "metient", "regularized parsimony", "normal parsimony", "constrained parsimony"]
structure=""

f1_violinplot(
    [mach2, metient, tlp_regularized, tlp_normal, constrained],
    fname_labels,
    structure,
)

f1_boxplot(
    [mach2, metient, tlp_regularized, tlp_normal, constrained],
    fname_labels,
    structure,
)
