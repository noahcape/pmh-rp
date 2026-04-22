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

def parse_multilabelings(fname, node_map):
    """
    Parse a file which contains labeling for vertices, return the labeling as well as the character_set
    Each line is a new labeling where the first word is the node label and the following words, seperated by tabs (or spaces) represent multi labels
    Return a dictionary where key = node and value = list of labels
    """
    labeling = {}
    character_set = []

    with open(fname, "r") as f:
        for line in f.readlines():
            split_line = line.split()

            if len(split_line) <= 1:
                raise ValueError(f"Lines in labeling file for clone tree ancestral reconstruction expected to be of the form <node> <label> (<label>)*")

            labels = split_line[1:]
            node = node_map[split_line[0]]
            labeling[node] = labels
            character_set.extend(labels)

    return (labeling, list(set(character_set)))
            

class MultiLabeledTree:
    def __init__(self, tree, labeling: dict, character_set):
        """
        Create a new multilabeled graph storing the graph and labels seperately
        
        Params:
            edgelist: graph edgelist as a list of tuples (source, dest)
            labeling: dictionary key = node, value = label(s); expect that label(s) are lists
        """        
        self.tree = tree
        self.labeling = labeling
        self.character_set = character_set

    @classmethod
    def from_tree(cls, tree, labeling: dict, character_set):
        """
        Create a new mutlilabeled tree from an already constructed tree
        """

        return cls(tree, labeling, character_set)
    
    @classmethod
    def from_edgelist(cls, edgelist, labeling: dict, character_set):
        tree = nx.DiGraph()
        tree.add_edges_from(edgelist)
        return cls(tree, labeling, character_set)

    def leaf_label_set(self):
        """
            Return set of labels appearing in leaves under v.
        """
        leaves = [u for u in self.tree.nodes if self.tree.out_degree(u) == 0]
        s = set()
        for u in leaves:
            labs = self.labeling.get(u, [])
            s.update(labs)
        return s

            

# labeling is a dict with node as entry and labels as the value
def format_metient_data(tree, labeling: dict, character_set, root_label, mutations=1):
    locations = character_set
    label_idx = {s: i for i, s in enumerate(locations)}

    rows = []

    for s in locations:
        for u in tree.nodes:
            rows.append(
                {
                    "anatomical_site_index": label_idx[s],
                    "anatomical_site_label": s,
                    "cluster_index": u,
                    "cluster_label": str(u).replace(":", "_"),
                    "present": (
                        1
                        if (u in labeling.keys() and s in labeling[u])
                        else 0
                    ),
                    "site_category": "primary" if s == root_label else "metastasis",
                    "num_mutations": mutations,
                }
            )

    return pd.DataFrame(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert tlp cancer evolution simulations into metient input"
    )
    parser.add_argument("tree", help="True tree in edgelist format")
    parser.add_argument("labeling", help="Leaf labeling of tree")
    parser.add_argument("-x", "--index", help="Index specifier for input", default=0)
    parser.add_argument("-o", "--output", help="Output file", default=0)
    parser.add_argument("-r", "--root", help="Root label", default=None)
    parser.add_argument("-m", "--metient-dir", help="Metient dir")
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

    if args.root == None:
        print("Must specify the root label for metient")    
        exit

    metient_dir = args.metient_dir
    try:
        os.mkdir(metient_dir)
    except Exception as _:
        print(f"Failed to make directory: {metient_dir}")
        exit

    tree = nx.read_edgelist(
        args.tree, nodetype=str, create_using=nx.DiGraph(), data=(("weight", float),)
    )

    node_map = {j: i for i, j in enumerate(tree.nodes)}
    print(node_map)
    edges = [(node_map[u], node_map[v]) for (u, v) in tree.edges]
    metient_tree = nx.from_edgelist(edges, create_using=nx.DiGraph())

    (labeling, character_set) = parse_multilabelings(args.labeling, node_map)
    multilabel_tree = MultiLabeledTree.from_tree(metient_tree, labeling, character_set)


    # Note that we pass metient tree which is the remapped nodes to 0-n indexing
    metient_df = format_metient_data(metient_tree, labeling,character_set, args.root)

    metient_data_path = f"{args.output}/metient.tsv"
    metient_tree_path = f"{args.output}/metient_tree.txt"

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
        with open(metient_tree_path, "w") as file:
            for u, v in metient_tree.edges:
                file.write(f"{u} {v}\n")
    except Exception as e:
        print(
            f"Failed to write tree edgelist as {args.tree} to {metient_tree_path} -- {e}."
        )
        exit

    print(f"Successfully wrote out metient tree format to {metient_tree_path}")

    print_config = metient.PrintConfig(visualize=False, verbose=False, k_best_trees=5)
    # set weight in decreasing order
    weights = metient.Weights(0.5, 0.4, 0.3)

    metient.evaluate_label_clone_tree(
        metient_tree_path,
        metient_data_path,
        weights,
        print_config,
        metient_dir,
        args.index,
        solve_polytomies=False,
        sample_size=10000,
    )
