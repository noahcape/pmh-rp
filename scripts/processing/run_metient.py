"""
Docstring for scripts.processing.create_metient_input


Take either mach2 simulations for clones trees or tlp simulations for leaf labeled as input
and process it into the input for metient. Since we only process one "patient" sample at a time
we only will have one character label - set to the seed used for the simulation.

Column name	            Description
-----------------------|--------------
anatomical_site_index	Zero-based index for anatomical_site_label column. Rows with the same anatomical site index and cluster_index will get pooled together.
anatomical_site_label	Name of the anatomical site
cluster_index	        If using a clustering method, the cluster index that this mutation belongs to. NOTE: this must correspond to the indices used in the tree txt file. Rows with the same anatomical site index and cluster_index will get pooled together.
cluster_label	        Name of the mutation or cluster of mutations. This is used in visualizations, so it should be short. NOTE: due to graphing dependencies, this string cannot contain colons.
present	                Must be one of 0 or 1. 1 indicates that this mutation/mutation cluster is present in this anatomical site, and 0 indicates that it is not.
site_category	        Must be one of primary or metastasis. If multiple primaries are specified, such that the primary label is used for multiple different anatomical site indices (i.e., the true primary is not known), we will run Metient multiple times with each primary used as the true primary. Output files are saved with the suffix _{anatomical_site_label} to indicate which primary was used in that run.
num_mutations	        The number of mutations in this cluster.


For number of mutations we should see what makes sense if not need to track this in the simulations
"""

import networkx as nx
import argparse
import pandas as pd
import os

from metient import metient


def format_metient_data(tree, labeling, root_label, mutations):
    locations = list(
        set([labeling.loc[leaf, "label"] for leaf in labeling.index])
    )
    nodes = [v for v in tree.nodes]
    leaves = [v for v in labeling.index]

    anatomical_site_idx = []
    anatomical_site_label = []
    cluster_idx = []
    cluster_label = []
    present = []
    site_category = []
    num_mutations = []
    label_idx = {i: j for j, i in enumerate(locations)}
    node_idx = {i: j for j, i in enumerate(nodes)}

    for s in label_idx:
        for u in node_idx:
            anatomical_site_idx.append(label_idx[s])
            anatomical_site_label.append(s)
            cluster_idx.append(node_idx[u])
            cluster_label.append(u)
            present.append(1 if int(u) in leaves and s == labeling.loc[int(u), "label"] else 0)
            site_category.append("primary" if s == root_label else "metastasis")
            num_mutations.append(5)


    metient_df = pd.DataFrame(
        {
            "anatomical_site_index": anatomical_site_idx,
            "anatomical_site_label": anatomical_site_label,
            "cluster_index": cluster_idx,
            "cluster_label": cluster_label,
            "present": present,
            "site_category": site_category,
            "num_mutations": num_mutations,
        }
    )

    return metient_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert tlp cancer evolution simulations into metient input"
    )
    parser.add_argument("tree", help="True tree in edgelist format")
    parser.add_argument("vertex_labeling", help="True vertex labeling in CSV format")
    parser.add_argument("leaf_labeling", help="Leaf labeling of tree")
    parser.add_argument("-o", "--output", help="output prefix")
    parser.add_argument("-x", "--index", help="Index specifier for input", default=0)
    return parser.parse_args()


if __name__ == "__main__":
    """
    Cluster idx = 0
    Cluster label = seed

    anatomical site index = node idx in tree
    anatomical site label = node's label
    present = 1
    site_category = primary iff it is the root in the labeling (thus need to parse the true labeling and specify - also then need to run the other methods with the true root information)
    num_mutations = n (check what this should be set as)
    """

    args = parse_args()

    parent = os.path.dirname(args.output)
    metient_dir = f"{parent}/metient"
    try:
        os.mkdir(metient_dir)
    except Exception as _:
        print(f"Failed to make directory: {metient_dir}")
        exit

    tree = nx.read_edgelist(
        args.tree, nodetype=str, create_using=nx.DiGraph(), data=(("weight", float),)
    )

    # for now lets assume that the weight is the number of mutations in this tree
    total_weight = sum(data["weight"] for _, _, data in tree.edges(data=True))

    try:
        full_labeling = pd.read_csv(args.vertex_labeling, sep=",").set_index("vertex")
    except Exception as _:
        column_names = ["vertex", "label"]
        full_labeling = pd.read_csv(
            args.vertex_labeling, sep="\t", header=None, names=column_names
        ).set_index("vertex")

    try:
        leaf_labeling = pd.read_csv(args.leaf_labeling, sep=",").set_index("leaf")
    except Exception as _:
        column_names = ["leaf", "label"]
        leaf_labeling = pd.read_csv(
            args.leaf_labeling, sep="\t", header=None, names=column_names
        ).set_index("leaf")

    num_leaves = len(leaf_labeling.index)
    root = [v for v in tree.nodes if tree.in_degree(v) == 0][0]
    root_label = full_labeling.loc[root, "label"]

    node_map = { j: i for i,j in enumerate(tree.nodes) }
    edges = [(node_map[u], node_map[v]) for (u,v) in tree.edges]
    metient_tree= nx.from_edgelist(edges, create_using=nx.DiGraph())

    metient_labeling = pd.DataFrame(
        {
            "leaf": [node_map[v] for v in leaf_labeling.index],
            "label": leaf_labeling["label"]
        }   
    ).set_index("leaf")

    metient_df = format_metient_data(metient_tree, metient_labeling, root_label, total_weight)

    metient_data_path = f"{args.output}_metient.tsv"
    metient_tree_path = f"{args.output}_metient_tree.txt"

    try:
        # write out the metient data to tsv
        metient_df.to_csv(metient_data_path, sep="\t", index=False)
    except Exception as e:
        print(
            f"Failed to write out metient data to {metient_data_path} with exception {e}"
        )
        exit

    print(f"Successfully wrote out metient data to {metient_data_path}")

    # convert the tree tsv file into a txt file - directly print the tree data to the file
    try:
        # tree_df = pd.read_csv(args.tree, sep="\t")
        # tree_df = tree_df.drop(tree_df.columns[2], axis=1)
        # tree_df.to_csv(metient_tree_path, sep=" ", index=False)
        with open(metient_tree_path, "w") as file:
            for (u,v) in metient_tree.edges:
                print(u,v)
                file.write(f"{u} {v}\n")
    except Exception as e:
        print(f"Failed to write tree edgelist as {args.tree} to {metient_tree_path} -- {e}.")
        exit

    print(f"Successfully wrote out metient tree format to {metient_tree_path}")

    print_config = metient.PrintConfig(visualize=True, verbose=False, k_best_trees=5)
    # set weight in decreasing order
    weights = metient.Weights(.5, .4, .3)

    metient.evaluate_label_clone_tree(
        metient_tree_path,
        metient_data_path,
        weights,
        print_config,
        metient_dir,
        args.index,
        solve_polytomies=False,
        sample_size=10000
    )
