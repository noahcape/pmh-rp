import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
import numpy as np
import os
from enum import Enum, auto

matplotlib.use("Agg")

matplotlib.rcParams.update({
    # General text
    "font.size": 14,
    "font.weight": "bold",

    # Axes labels and title
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "axes.titlesize": 18,
    "axes.titleweight": "bold",

    # Tick labels
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,

    # Legend
    "legend.fontsize": 13,
    "legend.title_fontsize": 14,

    # Grid
    "grid.linewidth": 1.2,
})

class FigureData(Enum):
    F1_SCORE = auto()
    RUNTIME = auto()
    PARSIMONY_GAP = auto()

DISPLAY_MAP = {
    "mach2": "MACH2",
    "metient": "Metient",
    "tlp_normal": "fastMACH",
    "tlp_polyclonal_tree": "fastMACH Tree",
    "tlp_polyclonal_dag": "fastMACH DAG",
    "tlp_tree": "fastMACH Tree",
    "tlp_dag": "fastMACH DAG",
    "tlp_reg": "LARCH",
}

def f1_score(df):
    """Compute F1 score from tp, fp, fn columns."""
    tp = df["tp"]
    fp = df["fp"]
    fn = df["fn"]
    denom = 2 * tp + fp + fn
    return (2 * tp) / denom.where(denom != 0)

def runtime(df):
    return df["time"]

def parsimony_gap(df):
    return df["ips"] - df["tps"]

def data_by_leaves_and_mean(dfs, labels, extractor, leaf_column="leaves", leaf_order=None):
    """
    Group data by leaf count, then sort methods within each leaf by mean F1 score.
    
    Returns:
        flat_scores: list of arrays (one per method per leaf group, non-empty)
        positions: x positions for each box/violin
        xticks: x-axis positions for leaf group labels
        xticklabels: leaf count labels
        method_colors: dict mapping method -> color
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Assign colors per method
    cmap = plt.get_cmap("tab10")
    method_colors = {m: cmap(i % 10) for i, m in enumerate(labels)}

    # Determine leaf order
    all_leaves = sorted(set(np.concatenate([df[leaf_column].unique() for (index, df) in enumerate(dfs) if leaf_column in dfs[index].columns])))
    if leaf_order:
        leaf_counts = [l for l in leaf_order if l in all_leaves]
    else:
        leaf_counts = all_leaves

    flat_scores = []
    positions = []
    xticks = []
    xticklabels = []

    pos = 1  # starting x position

    for leaf in leaf_counts:
        # Gather method scores for this leaf count
        method_scores = []
        for df, label in zip(dfs, labels):
            if not df.empty:
                scores = extractor(df[df[leaf_column] == leaf]).dropna().values
                if len(scores) != 0:
                    method_scores.append((label, scores))  # keep only non-empty arrays

        if not method_scores:
            continue  # skip leaf group if no method has data

        # Add scores and positions
        start_pos = pos
        for label, scores in method_scores:
            flat_scores.append((label, scores))
            positions.append(pos)
            pos += 1
        end_pos = pos - 1

        # Center x-axis label for this leaf group
        group_center = (start_pos + end_pos) / 2
        xticks.append(group_center)
        xticklabels.append(str(leaf))

        pos += 1  # extra space between leaf groups

    return flat_scores, positions, xticks, xticklabels, method_colors

def boxplot(dfs, labels, simulation_type, extractor, title, data_type: FigureData, out: Path, leaf_column="leaves", leaf_order=None):
    flat_scores, positions, xticks, xticklabels, method_colors = data_by_leaves_and_mean(
        dfs, labels, extractor, leaf_column, leaf_order
    )

    scores_only = [scores for label, scores in flat_scores]

    min_pos = min(positions)
    positions = [(p - min_pos) * 0.15 for p in positions]
    xticks = [(x - min_pos) * 0.15 for x in xticks]
    plt.figure(figsize=(max(8, len(flat_scores) * 0.5), 5))
    box = plt.boxplot(
        scores_only,
        positions=positions,
        widths=0.10,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=2),
    )
    plt.margins(x=0)
    # plt.boxplot(data, widths=0.4)

    # Color boxes and whiskers/caps by actual method label
    for patch, (label, _) in zip(box["boxes"], flat_scores):
        patch.set_facecolor(method_colors[label])
        patch.set_alpha(0.7)

    for i, (label, _) in enumerate(flat_scores):
        color = method_colors[label]
        box["whiskers"][2*i].set_color(color)
        box["whiskers"][2*i + 1].set_color(color)
        box["whiskers"][2*i].set_linewidth(1.5)
        box["whiskers"][2*i + 1].set_linewidth(1.5)
        box["caps"][2*i].set_color(color)
        box["caps"][2*i + 1].set_color(color)
        box["caps"][2*i].set_linewidth(1.5)
        box["caps"][2*i + 1].set_linewidth(1.5)
        flier = box["fliers"][i]
        flier.set_markerfacecolor(color)
        flier.set_markeredgecolor(color)
        flier.set_alpha(0.7)
        flier.set_markersize(4)

    match data_type:
        case FigureData.RUNTIME:
            plt.yscale("log")
            plt.ylabel(f"{title} - log(sec)")
        case FigureData.F1_SCORE:
            plt.ylabel(title)
            plt.ylim(-0.05, 1.05)
        case FigureData.PARSIMONY_GAP:
            plt.ylabel(f"{title} (inferred - true)")

    # plt.title(f"{title} Distribution Comparison ")
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.xticks(ticks=xticks, labels=xticklabels)
    plt.xlabel("Leaves")
    plt.xlim(min(positions) - 0.2, max(positions) + 0.2)

    # Legend
    for method, color in method_colors.items():
        plt.plot([], [], color=color, linewidth=8, label=DISPLAY_MAP[method])

    match data_type:
        case FigureData.F1_SCORE:
            plt.legend(title="Method", loc="lower right")
        case FigureData.PARSIMONY_GAP:
            plt.ylabel(title)
    
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("Saved:", out.resolve())
    plt.close()


def violinplot(dfs, labels, simulation_type, extractor, title, runtime=False, leaf_column="leaves", leaf_order=None):
    flat_scores, positions, xticks, xticklabels, method_colors = data_by_leaves_and_mean(
        dfs, labels, extractor, leaf_column, leaf_order
    )

    plt.figure(figsize=(max(8, len(flat_scores) * 0.5), 5))

    for xi, (label, scores) in zip(positions, flat_scores):
        color = method_colors[label]
        if len(scores) == 0:
            continue

        parts = plt.violinplot([scores], positions=[xi], showmeans=True, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        for key in ["cbars", "cmins", "cmaxes", "cmedians", "cmeans"]:
            parts[key].set_color("black")
            parts[key].set_linewidth(1.5)

    if runtime:
        plt.yscale("log")
        plt.ylabel(f"{title} - log(sec)")
    else:
        plt.ylabel(title)
        plt.ylim(-0.05, 1.05)

    plt.title(f"{title} Distribution Comparison ({simulation_type})")
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.xticks(ticks=xticks, labels=xticklabels)
    plt.xlabel("Number of Leaves")

    # Legend
    for method, color in method_colors.items():
        plt.plot([], [], color=color, linewidth=8, label=method)
    
    if not runtime:
        plt.legend(title="Method", loc="lower right")

    plt.tight_layout()
    underscored_title = title.lower().replace(" ", "_")
    plt.savefig(f"{underscored_title}_violin_{simulation_type}.png", dpi=300, bbox_inches="tight")
    print("Saved:", Path(f"{underscored_title}_violin_{simulation_type}.png").resolve())
    plt.close()
    
def visualize_simulation_results(path: Path, generations: list[int]):
    results = defaultdict(list)
    def load_csv(p: Path):
        if os.path.exists(p):
            return pd.read_csv(p).dropna(how="all")
        else: return None
        
    def score_file(method: str): return f"{method}_results.csv"
    
    mach2 = "mach2"
    metient = "metient"
    tlp_normal = "tlp_normal"
    tlp_reg = "tlp_reg"
    tlp_dag = "tlp_dag"
    tlp_tree = "tlp_tree"
    
    methods = [mach2, metient, tlp_normal, tlp_dag, tlp_tree, tlp_reg]
    methods = [mach2, metient, tlp_normal, tlp_reg]
    
    for g in generations:
        base_dir = path / str(g)
        print(g)
        for m in methods:
            score_file_name = base_dir / score_file(m)
            df = load_csv(score_file_name)
            
            if df is not None:
                df["leaves"] = np.pow(2, g)
                results[m].append(df)
                
    concatenated_results = [
            pd.concat(results[method], ignore_index=True) 
            if len(results[method]) != 0 else pd.DataFrame() 
            for method in methods
        ]

    boxplot(
        concatenated_results,
        methods,
        "",
        f1_score,
        "F1 Score",
        FigureData.F1_SCORE,
        Path("figures/sims_results") / "rust_sims_v2_f1_score.png"        
    )

    boxplot(
        concatenated_results,
        methods,
        "",
        runtime,
        "Runtime",
        FigureData.RUNTIME,
        Path("figures/sims_results") / "rust_sims_v2_runtime.png"
    )
    
    return None

def visualize_error_simulations():
    migration_rate = [0.002]
    # leaves = [50, 75, 100]
    leaves = [50, 75, 100, 500, 750, 1000]
    flips = [5, 10]
    constraints = ["polyclonal_dag"]
    # constraints = ["none", "polyclonal_dag", "polyclonal_tree"]

    base = Path("/n/fs/ragr-research/projects/pmh-rp/sims")

    dfs = {
        "normal": defaultdict(list),
        "perturbed": defaultdict(list),
        "flip=10": defaultdict(list),
        "flip=5": defaultdict(list)
    }


    def load_csv(path):
        if os.path.exists(path):
            return pd.read_csv(path).dropna(how="all")
        else:
            return None
        
    mach2 = "mach2"
    metient = "metient"
    tlp_normal = "tlp_normal"
    tlp_reg = "tlp_reg"
    tlp_dag = "tlp_polyclonal_dag"
    tlp_tree = "tlp_polyclonal_tree"
        
    
    # go in this order same as the CP_ clone data
    # methods = ["tlp_normal", "tlp_dag", "tlp_tree", "mach2", "metient", "regularized"]
    # methods = [mach2, metient, tlp_normal, tlp_dag, tlp_reg]
    methods = [mach2, metient, tlp_normal, tlp_dag, tlp_reg]
    

    def score_file(method: str, error_style):
        if method == "mach2" or method == "metient":
            return f"{method}_{error_style}.csv"
        if method.startswith("tlp"):
            if error_style == "normal":
                return f"{method}.csv"
            else:
                return f"{method}_{error_style}.csv"
    
    for mrate in migration_rate:
        for num_leaves in leaves:
            for flip in flips:
                for constraint in constraints:
                    params = lambda i : f"leaves={num_leaves}_constraint={constraint}_mrate={mrate}_flip={i}"
                    root = base / params(flip)

                    update_dfs = ["normal", "perturbed", f"flip={flip}"]

                    for method in methods:
                        for df_name in update_dfs:
                            score_file_name = root / score_file(method, df_name)
                            df = load_csv(score_file_name)

                            if df is not None:
                                df["leaves"] = num_leaves
                                dfs[df_name][method].append(df)

    concatenated_dfs = {
        # "normal": [pd.concat(dfs["normal"][method], ignore_index=True) if len(dfs["normal"][method]) != 0 else pd.DataFrame() for method in methods],
        "perturbed": [pd.concat(dfs["perturbed"][method], ignore_index=True) if len(dfs["perturbed"][method]) != 0 else pd.DataFrame() for method in methods],
        # "flip=10": [pd.concat(dfs["flip=10"][method], ignore_index=True) if len(dfs["flip=10"][method]) != 0 else pd.DataFrame() for method in methods],
        # "flip=5": [pd.concat(dfs["flip=5"][method], ignore_index=True) if len(dfs["flip=5"][method]) != 0 else pd.DataFrame() for method in methods],
    }

    # maybe merge dag and tree into a constrained version
    for key in concatenated_dfs:
        title = f""
        boxplot(
            concatenated_dfs[key],
            methods,
            title,
            f1_score,
            "F1 Score",
            FigureData.F1_SCORE
        )

        boxplot(
            concatenated_dfs[key],
            methods,
            title,
            runtime,
            "Runtime",
            FigureData.RUNTIME
        )
        # boxplot(
        #     concatenated_dfs[key],
        #     methods,
        #     key,
        #     parsimony_gap,
        #     "Parsimony Gap",
        #     FigureData.PARSIMONY_GAP
        # )

    return None

def visualize_normal_simulations():
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

                # mach2
                df = load_csv(root / f"mach2_{constraint}.csv")
                df["leaves"] = num_leaves
                dfs["mach2"].append(df)

                # metient (only for small leaf counts)
                if num_leaves < 100:
                    df = load_csv(root / "metient.csv")
                    df["leaves"] = num_leaves
                    dfs["metient"].append(df)

                # tlp_regularized
                df = load_csv(root / f"tlp_{constraint}_r=1.csv")
                df["leaves"] = num_leaves
                dfs["tlp_regularized"].append(df)

                # tlp_normal
                if constraint == "none":
                    df = load_csv(root / "tlp_none_r=0.csv")
                else:
                    df = load_csv(root / f"tlp_{constraint}.csv")
                df["leaves"] = num_leaves
                dfs["tlp_normal"].append(df)

                # constrained tlp
                if constraint == "polyclonal_tree":
                    df = load_csv(root / f"tlp_{constraint}_r=0.csv")
                    df["leaves"] = num_leaves
                    dfs["tlp_tree"].append(df)
                elif constraint == "polyclonal_dag":
                    df = load_csv(root / f"tlp_{constraint}_r=0.csv")
                    df["leaves"] = num_leaves
                    dfs["tlp_dag"].append(df)


    # Final concatenation
    mach2 = pd.concat(dfs["mach2"], ignore_index=True)
    metient = pd.concat(dfs["metient"], ignore_index=True)
    tlp_regularized = pd.concat(dfs["tlp_regularized"], ignore_index=True)
    tlp_normal = pd.concat(dfs["tlp_normal"], ignore_index=True)
    tlp_tree = pd.concat(dfs["tlp_tree"], ignore_index=True)
    tlp_dag = pd.concat(dfs["tlp_dag"], ignore_index=True)
    constrained = pd.concat([tlp_tree, tlp_dag], ignore_index=True)

    fname_labels = ["metient", "mach2", "normal parsimony", "constrained parsimony", "regularized parsimony"]
    structure=""

    # f1_violinplot(
    #     [mach2, metient, tlp_regularized, tlp_normal, constrained],
    #     fname_labels,
    #     structure,
    #     f1_score,
    #     "Migration graph f1 score"
    # )

    boxplot(
        [metient, mach2, tlp_normal, constrained, tlp_regularized],
        fname_labels,
        structure,
        f1_score,
        "Migration graph f1 score",
        FigureData.F1_SCORE
    )

    boxplot(
        [metient, mach2, tlp_normal, constrained, tlp_regularized],
        fname_labels,
        structure,
        runtime,
        "Runtime",
        FigureData.RUNTIME
    )

    # boxplot(
    #     [tlp_normal, constrained],
    #     ["normal parsimony", "constrained parsimony"],
    #     structure,
    #     runtime,
    #     "Runtime test",
    #     True
    # )


if __name__ == "__main__":
    base = Path("/n/fs/ragr-research/projects/pmh-rp")
    visualize_simulation_results(base / "rust_sims_v2", [6,7,8,9,10])