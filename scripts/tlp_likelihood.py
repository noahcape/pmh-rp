import sys
import argparse
import itertools
import json

import pandas as pd
import numpy as np
import networkx as nx
import pyomo.environ as pyo
from scipy.linalg import expm

from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict
from gurobipy import GRB


"""Find the transition cost along edge u,v between sites a,b"""


def edge_cost(e, a, b, Qt_map):
    if a is None or b is None:
        return 0

    if a == b:
        return 0
    else:
        return Qt_map[e][a][b]


"""Precompute all edge dependent weights"""


def precompute_Qt(T: nx.DiGraph, Q):
    Qt_map = {}
    for u, v, m in T.edges(data=True):
        t = m["weight"]
        print((u, v))
        try:
            Qt_map[(u, v)] = expm(Q * t)
        except Exception as _:
            raise ValueError("Tree must have branch lengths")

    return Qt_map


"""Validate that Q is a rate matrix"""


def validate_Q(Q: np.ndarray):
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
        raise ValueError("Rate matrix must be a square matrix.")

    if np.any(Q.sum(axis=1)):
        raise ValueError("Rate matrix rows must sum to zero.")


def process_model(model, tree, character_set, vertex_labeling={}):
    if len(vertex_labeling.keys()) == 0:
        # compute (an) optimal vertex label
        vertex_labeling = {}
        for u, v in tree.edges:
            for c1, c2 in itertools.product(character_set, character_set):
                if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                    vertex_labeling[u] = c1
                    vertex_labeling[v] = c2

    # compute the new migration multi-graph
    migration_graph = defaultdict(int)
    for u, v in tree.edges:
        if u not in vertex_labeling or v not in vertex_labeling:
            continue

        if vertex_labeling[u] != vertex_labeling[v]:
            migration_graph[(vertex_labeling[u], vertex_labeling[v])] += 1

    return migration_graph, vertex_labeling


def write_results(migration_graph, vertex_labeling, out, params={}):
    logger.info("Writing out results from model solution.")

    migrations = 0
    migation_pattern_edges = 0

    # writes a new optimal labeling to a file
    with open(f"{out}_vertex_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for node, label in vertex_labeling.items():
            f.write(f"{node},{label}\n")

    # write a new migration multi-graph to a file
    with open(f"{out}_migration_graph.csv", "w") as f:
        f.write("source,target,count\n")
        for (i, j), count in migration_graph.items():
            migrations += count
            migation_pattern_edges += 1

            f.write(f"{i},{j},{count}\n")

    # write the objective value to a file (json)
    with open(f"{out}_results.json", "w") as f:
        results = {}

        results["migrations"] = migrations
        results["migration_pattern_edges"] = migation_pattern_edges
        # insert parameters of the lp
        results.update(params)
        f.write(json.dumps(results))

    print(
        f"<------ Migrations in labeling: {migrations} ------>\n<------ Number of edges in migration graph: {migation_pattern_edges} ------>"
    )


def set_objective(
    T, character_set, model, dist_f, regularization_param, add_migrations=True
):
    """
    Set the objective of the model to be the sum of the number of edges in the migration pattern
    and the number of migrations in the tree.
    """
    logger.info(f"Setting object with weight: {regularization_param}")

    if add_migrations:
        # append migrations since we need it in our objective
        append_migrations(model, tree, character_set)

    # set objective num_edges + num_migrations
    num_edges = sum(
        model.migrations[(i, j)] for i in character_set for j in character_set if i != j
    )

    # \sum_{u,v} \sum_{c,c'} x_{u,v,c,c'} * d((u,v), c,c')
    migrations = sum(
        model.decisions[u, v, c, c2] * dist_f((u, v), c, c2)
        for u, v in T.edges
        for c in character_set
        for c2 in character_set
    )

    model.objective_expr = (regularization_param * num_edges) + migrations
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
    return model


"""
Creates the tree labeling polytope whose vertices correspond
to the set of solutions to the (unconstrained) maximum parsimony 
problem.
"""


def create_tree_labeling_polytope(
    T, root, character_set, label_f, dist_f, root_label=None
):
    model = pyo.ConcreteModel()

    logger.info(
        f"Creating tree labeling polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges."
    )
    logger.info(f"Character set has size: {len(character_set)}")

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()
    # set domain to be [0, 1]

    # model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.Binary, initialize=1)
    model.decisions = pyo.Var(
        T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0
    )

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w,c
    model.edge_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                model.edge_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set)
                    - sum(model.decisions[v, w, c, c2] for c2 in character_set)
                    == 0
                )

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), v is a leaf, and c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), v is a leaf, and c' is not the label of v
    model.leaf_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        if label_f(v) == None:
            continue
        for c in character_set:
            if c == label_f(v):
                model.leaf_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set) == 1
                )
            else:
                model.leaf_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set) == 0
                )

    # set root label if provided
    if root_label is not None:
        model.root_constraint = pyo.ConstraintList()
        model.root_constraint.add(
            sum(
                model.decisions["dummy_root", root, c, root_label]
                for c in character_set
            )
            == 1
        )

    return model


"""
Appends migration variables to the tree labeling polytope.
"""


def append_migrations(model, T, character_set):
    logger.info("Adding migration constraints.")
    edges = [(i, j) for i in character_set for j in character_set if i != j]
    model.migrations = pyo.Var(edges, domain=pyo.Binary, initialize=1)
    model.migration_constraints = pyo.ConstraintList()
    for u, v, c1, c2 in model.decisions:
        if c1 == c2:
            continue
        model.migration_constraints.add(
            model.decisions[u, v, c1, c2] <= model.migrations[c1, c2]
        )


"""
Solves a generalization of the MACHINA parsimonious migration 
history problem using the tree labeling polytope.
"""


def fast_machina(tree, character_set, leaf_f, dist_f, root, args, mip_gap=0.15):
    tree = tree.copy()

    """ Step 1: Setup the MILP using Gurobi callbacks, if necessary """
    solver = pyo.SolverFactory("gurobi_persistent")

    model = create_tree_labeling_polytope(
        tree, root, character_set, leaf_f, dist_f, root_label=args.label
    )
    # set objective
    model = set_objective(tree, character_set, model, dist_f, args.regularizer_weight)
    solver.set_instance(model)

    """ Step 3: Solve the full MILP """
    logger.info("Solving full MILP model")
    solver.options["MIPGap"] = mip_gap
    solver.solve(model, tee=True, warmstart=True)

    # compute (an) optimal vertex labeling
    vertex_labeling = {}
    for u, v in tree.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2

    return vertex_labeling, model.objective()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Constrained tree labeling using the tree labeling polytope."
    )

    parser.add_argument("tree", help="Tree in edgelist format")
    parser.add_argument("labels", help="Leaf labeling as a CSV file")
    parser.add_argument("qmatrix", help="Q rate matrix in csv format")

    parser.add_argument("-o", "--output", help="Output prefix", default="result")
    parser.add_argument("-r", "--root", help="Root of the tree", default=None)
    parser.add_argument("-l", "--label", help="Root label", default=None)
    parser.add_argument(
        "-w",
        "--regularizer-weight",
        help="likelihood regularization weight",
        type=float,
        default=1.0,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(
        args.tree, create_using=nx.DiGraph(), data=(("weight", float),)
    )

    # fall back to normal parsing for phylogenetic (non-clone tree) tree instance
    labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")

    # load the character set as the set of all leaf labels
    character_set = labels_csv[labels_csv["label"] != "None"]["label"].unique()
    character_set_idx = {label: i for i, label in enumerate(character_set)}

    if not nx.is_directed_acyclic_graph(tree):
        raise ValueError("Graph is not a tree, it contains cycles.")

    if not nx.is_weakly_connected(tree):
        raise ValueError("Graph is not connected, it is a forest.")

    # process root label
    if not hasattr(args, "label"):
        args.label = None

    if args.label is not None and args.label not in character_set:
        logger.warning(f"Root label {args.label} not in character set, removing it.")
        args.label = None

    if args.root is not None:
        root = args.root
    else:
        roots = [node for node in tree.nodes if len(list(tree.predecessors(node))) == 0]
        if len(roots) != 1:
            raise ValueError(f"Tree has {len(roots)} roots, please specify the root.")
        root = roots[0]

    if "dummy_root" in tree.nodes:
        raise ValueError("Cannot use dummy_root in nodes")

    # add a dummy root node
    tree.add_node("dummy_root")
    tree.add_edge("dummy_root", root, weight=0.0)

    Q = np.loadtxt(args.qmatrix, delimiter=",")
    # validate that Q is a rate matrix
    validate_Q(Q)
    # precompute edge dependent costs
    Qt_map = precompute_Qt(tree, Q)

    # defines the distance function between characters x and y along an edge e
    def dist_f(e, x, y):
        a = character_set_idx[x]
        b = character_set_idx[y]
        return edge_cost(e, a, b, Qt_map)

    def is_int(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    def is_leaf(T, node):
        return len(T[node]) == 0

    # defines the leaf labeling function
    def leaf_f(node):
        if is_int(node):
            node = int(node)

        if node not in labels_csv.index:
            # print("Not found", node)
            return None

        y = labels_csv.loc[node, "label"]
        if y == "None":
            return None

        return y

    if not all(
        [leaf_f(node) is not None for node in tree.nodes if is_leaf(tree, node)]
    ):
        unlabeled_leaves = [
            node for node in tree.nodes if is_leaf(tree, node) and leaf_f(node) is None
        ]
        raise ValueError(f"Leaves {unlabeled_leaves} are unlabeled.")

    if len(character_set) == 1:
        logger.warning("Character set has size 1, inferring trivial labeling.")
        vertex_labeling = {node: character_set[0] for node in tree.nodes}
        lp_obj = 0
        obj = 0
    else:
        # computes the vertex labeling using the specified method
        lp_obj = None
        vertex_labeling, obj = fast_machina(
            tree, character_set, leaf_f, dist_f, root, args
        )

    (migration_graph, vertex_labeling) = process_model(
        None, tree, character_set, vertex_labeling
    )
    write_results(migration_graph, vertex_labeling, args.output)
