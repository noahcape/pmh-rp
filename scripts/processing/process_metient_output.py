"""
    read and process metient output from *.pkl.gz

    Pkl key name	        Description
    anatomical_sites	    a list of anatomical sites in the order used for the matrices detailed below.
    node_info	            list of dictionaries, in order from best to worst solution. This is solution specific because reolving polytomies can change the tree. Each dictionary maps node index (as used for the matrices detailed below) to a tuple: (label, is_leaf, is_polytomy_resolver_node) used on the tree. The reason labels can be different from what is inputted into Metient is that Metient adds leaf nodes which correspond to the inferred presence of each node in anatomical sites. Each leaf node is labeled as <parent_node_name>_<anatomical_site>.
    node_labels	            list of numpy ndarrays, in order from best to worst solution. Each numpy array is a matrix (shape: len(ordered_anatomical_sites), len(node_info[x])), where x is the xth best solution. Row i corresponds to the site at index i in ordered_anatomical_sites, and column j corresponds to the node with label node_info[x][j][0]. Each column is a one-hot vector representing the location inferred by Metient for that node.
    parents	                list of numpy 1-D arrays, in order from best to worst tree. Each is a an array (shape: len(node_info[x])), where x is the xth best solution. The value at index i is the parent of node i. The root node will have a -1 at its index.
    observed_proportions	numpy ndarray (shape: len(ordered_anatomical_sites), num_clusters). Row i corresponds to the site at index i in ordered_anatomical_sites, and column j corresponds to the node with label node_info[x][j][0]. A value at i,j greater than 0.05 indicates that that node is present in that antomical site. These are the nodes that get added as leaf nodes.
    losses	                a list of the losses, from best to worst solution.
    probabilities	        a list of the probabilities, from best to worst solution.
    primary_site	        str, the name of the anatomical site used as the primary site.
    loss_info	            a list of the dicts, from best to worst solution. Each dictionary contains the unweighted components of the loss (e.g. migration number, comigration number, etc.)

"""

import os
import gzip
import pickle
import subprocess
import argparse

def full_node_map(data, solution_idx=0):
    anatomical_sites = data["anatomical_sites"]
    node_info = data["node_info"][solution_idx]
    node_labels = data["node_labels"][solution_idx]
    parents = data["parents"][solution_idx]

    result = {}

    for idx, (label, is_leaf, is_poly) in node_info.items():
        site = anatomical_sites[node_labels[:, idx].argmax()]
        parent = parents[idx] if parents[idx] != -1 else None

        result[idx] = {
            "idx": idx,
            "label": label,
            "site": site,
            "parent": parent,
            "is_leaf": is_leaf,
            "is_polytomy_resolver": is_poly,
        }

    return result

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process and score metient output"
    )
    parser.add_argument("metient", help="base directory to find files") # this is going to be the metient dir
    parser.add_argument("timing", help="Timing file from running metient")
    parser.add_argument("tree", help="Ground truth tree")
    parser.add_argument("solved_tree", help="Tree solved reconstruction on")
    parser.add_argument("labeling", help="Full labeling of tree")
    parser.add_argument("-s", "--stats", help="File to collect stats", default=None, required=False)
    return parser.parse_args()


if __name__ == "__main__":
    score_labeling_file = "/n/fs/ragr-research/projects/pmh-rp/scripts/processing/score_result.py"
    draw_labeling_file = "/n/fs/ragr-research/projects/pmh-rp/scripts/plots/draw_colored_tree.py"
    extract_stats = "/n/fs/ragr-research/projects/pmh-rp/extract_stats.sh"
    
    args = parse_args()

    metient_dir = args.metient
    timing_file = args.timing
    tree = args.tree
    solved_on_tree = args.solved_tree
    full_labeling = args.labeling

    # put the score file inside the metient dir
    score_file = lambda i : f"{metient_dir}/{i}_metient_results.json"
    # put the labeling file which will be created inside the metient dir
    labeling_file = lambda i : f"{metient_dir}/{i}_metient_labeling.csv"
    tree_file = lambda i : f"{metient_dir}/{i}_metient_edgelist.tsv"
    
    # get the metient output
    pkl_files = [
        f for f in os.listdir(metient_dir)
        if os.path.isfile(os.path.join(metient_dir, f)) and f.endswith(".pkl.gz")
    ]

    if len(pkl_files) == 0:
        print("No results found.")
        exit(1)

    metient_output = pkl_files[0]

    with gzip.open(f"{metient_dir}/{metient_output}", "rb") as f:
        
        data = pickle.load(f)
        
        num_solutions = len(data["node_info"])
        
        for i in range(num_solutions):
            edge_list = []
            label_map = []
            node_map = full_node_map(data, i)
            
            for node in node_map:
                entry = node_map[node]

                label = entry["label"]
                site = entry["site"]
                parent = entry["parent"]
                is_leaf = entry["is_leaf"]
                idx = entry["idx"]
            
                node_idx = idx
                parent_idx = int(parent) if not parent == None else None
                label = site
                
                if not parent_idx == None and parent_idx != node_idx:
                    edge_list.append((parent_idx, node_idx))
                    
                if (node_idx, label) not in label_map:
                    label_map.append((node_idx, label))
                
            # labeling csv with vertex,label (s<vertex>)
            with open(labeling_file(i), "w") as file:
                file.write("vertex,label\n")
                for (node_idx, label) in label_map:
                    file.write(f"{node_idx},{label}\n")

            with open(tree_file(i), "w") as file:
                for (u, v) in edge_list:
                    file.write(f"{u}\t{v}\n")

            # subprocess.run(
            #         [
            #             "conda", "run", "-n", "tlp",
            #             "python",
            #             draw_labeling_file,
            #             tree_file(i),
            #             labeling_file(i),
            #             "-o", f"{metient_dir}/metient_{i}",
            #             "-s"
            #         ],
            #         check=True
            #     )

            # subprocess.run(
            #     [
            #         "dot", "-Tpng", "-Gdpi=300",
            #         f"{metient_dir}/metient_{i}_color_graph.dot", "-o", f"{metient_dir}/metient_{i}_color_graph.png"
            #     ],
            #     check=True
            # )
            # subprocess.run(
            #     [
            #         "dot", "-Tpng", "-Gdpi=300",
            #         f"{metient_dir}/metient_{i}_colored_tree.dot", "-o", f"{metient_dir}/metient_{i}_color_tree.png"
            #     ],
            #     check=True
            # )

            if args.stats != None:
                # score result
                subprocess.run(
                    [
                        "conda", "run", "-n", "tlp",
                        "python",
                        score_labeling_file,
                        solved_on_tree,
                        tree,
                        full_labeling,
                        labeling_file(i),
                        timing_file,
                        "-o", score_file(i),
                    ],
                    check=True
                )
                
                #extract stats
                subprocess.run(
                    [
                        "bash", extract_stats, score_file(i), args.stats
                    ],
                    check=True
                )
