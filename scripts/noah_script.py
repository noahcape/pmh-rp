import pandas as pd
import matplotlib.pyplot as plt


def jaccard_index_boxplot(dfs, labels, simulation_type):
    indices = [df["ji"] for df in dfs]

    plt.figure(figsize=(7, 5))
    plt.boxplot(
        indices,
        tick_labels=labels,
    )

    plt.ylabel("Jaccard Index")
    plt.title(f"Jaccard Index Distribution Comparison ({simulation_type})")
    plt.ylim(top=1.1, bottom=-0.1)
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()


def parsimony_score_gap(fnames, labels, simulation_type):
    gaps = [
        ((pd.read_csv(f)["tps"] - pd.read_csv(f)["ips"]).abs() * 100)
        / pd.read_csv(f)["tps"]
        for f in fnames
    ]
    print(gaps)
    plt.figure(figsize=(7, 5))
    plt.boxplot(
        gaps,
        tick_labels=labels,
    )
    plt.ylabel("Percentage Score Gap (true & inferred)")
    plt.title(f"Parsimony Score Gap ({simulation_type})")
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()


def f1_boxplot(dfs, labels, simulation_type):
    def f1_score(df):
        fp = df["fp"]
        fn = df["fn"]
        tp = df["tp"]

        f1 = (2 * tp) / (2 * tp + fp + fn)

        return f1

    f1_scores = [f1_score(df) for df in dfs]

    plt.figure(figsize=(7, 5))
    plt.boxplot(f1_scores, tick_labels=labels)

    plt.ylabel("F1 Score")
    plt.ylim(top=1.1, bottom=-0.1)
    plt.title(f"F1 Score Distribution Comparison ({simulation_type})")
    plt.grid(axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()


# tlp_stats = ["./data/clone_trees/sims/m5/M/tlp_stats.csv", "./data/clone_trees/sims/m5/R/tlp_stats.csv", "./data/clone_trees/sims/m5/S/tlp_stats.csv"]
# mach2_stats = ["./data/clone_trees/sims/m5/M/mach2_stats.csv", "./data/clone_trees/sims/m5/R/mach2_stats.csv", "./data/clone_trees/sims/m5/S/mach2_stats.csv"]
tlp_normal_stats = "./examples/cancer_evolution/sims_tree_large/tlp_normal_stats.csv"
tlp_tree_stats = "./examples/cancer_evolution/sims_tree_large/tlp_tree_stats.csv"
tlp_stats = "./examples/cancer_evolution/sims_tree_large/tlp_stats.csv"
mach2_stats = "./examples/cancer_evolution/sims_tree_large/mach2_stats.csv"

# tlp_dfs = pd.concat([pd.read_csv(f) for f in tlp_stats], ignore_index=True)
# mach2_dfs = pd.concat([pd.read_csv(f) for f in mach2_stats], ignore_index=True)
tlp_normal_dfs = pd.read_csv(tlp_normal_stats)
tlp_tree_dfs = pd.read_csv(tlp_tree_stats)
tlp_dfs = pd.read_csv(tlp_stats)
mach2_dfs = pd.read_csv(mach2_stats)

fname_labels = ["normal parsimony", "tree constrained parsimony","regularized parsimony", "mach2"]
structure="mig_rate=3e-3 leaves=50"

f1_boxplot(
    [tlp_normal_dfs, tlp_tree_dfs, tlp_dfs, mach2_dfs],
    fname_labels,
    structure,
)
# for structure in ["polyclonal_tree", "polyclonal_dag"]:
#     fname1 = f"./examples/cancer_evolution/simulations/simulations_{structure}_results_large_no_reg.csv"
#     fname2 = f"./examples/cancer_evolution/simulations/simulations_{structure}_results_large_regularized.csv"
#     fname3 = f"./examples/cancer_evolution/simulations/simulations_{structure}_results_{structure}.csv"

#     # fnames = [fname1, fname2, fname3]
#     # fname_labels = ["Normal Parsimony", "Regularized Parsimony", "Constrained Parsimony"]
#     fnames = [fname1, fname3]
#     fname_labels = ["Normal Parsimony", "Constrained Parsimony"]


#     parsimony_score_gap(
#         fnames,
#         fname_labels,
#         structure,
#     )
    # jaccard_index_boxplot(
    #     fnames,
    #     fname_labels,
    #     structure,
    # )
    # f1_boxplot(
    #     fnames,
    #     fname_labels,
    #     structure,
    # )
