import argparse
import itertools
import json
import sys

import pandas as pd
import numpy as np
import networkx as nx

from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict

import pyomo.environ as pyo
from gurobipy import GRB
import gurobipy as gp


def parse_multilabelings(fname):
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
                raise ValueError(
                    f"Lines in labeling file for clone tree ancestral reconstruction expected to be of the form <node> <label> (<label>)*"
                )

            labels = split_line[1:]
            labeling[split_line[0]] = labels
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

    def augment_tree(self, out=None) -> tuple[nx.DiGraph, dict, list, list]:
        """
        Augment the tree into a DAG where multi-labeled vertices are expanded
        into a 'search graph' encoding possible resolutions
        """
        augmented_labeling = {}
        restricted_edges = []
        search_complexes = []
        N: nx.DiGraph = self.tree.copy()

        for node in self.tree.nodes:
            labels = self.labeling.get(node, [])

            if len(labels) > 1:
                parents = N.predecessors(node)
                children = N.successors(node)
                new_nodes = [f"{node}_{i}" for i in range(len(labels))]
                new_edges = [
                    (new_nodes[i], new_nodes[j])
                    for i in range(len(labels))
                    for j in range(i + 1, len(labels))
                ]
                in_edges = [(parent, new_nodes[0]) for parent in parents]
                # iterate over children first becuase it is an iterator and if if its second it gets consumed on the first iteration
                out_edges = [
                    (new_nodes[i], child)
                    for child in children
                    for i in range(len(labels))
                ]
                N.add_edges_from(new_edges)
                N.add_edges_from(in_edges)
                N.add_edges_from(out_edges)
                N.remove_node(node)

                search_complex = (new_nodes, labels)
                search_complexes.append(search_complex)
                restricted_edges.extend(new_edges)

                # add these nodes to the augmented labeing dictionary for later constraint
                for new_node in new_nodes:
                    augmented_labeling[new_node] = labels

        if out is not None:
            with open(f"{out}_augmented_network.edgelist", "w+") as f:
                for u, v in N.edges:
                    f.write(f"{u}\t{v}\n")

        return (N, augmented_labeling, restricted_edges, search_complexes)

    def dist_f(self, a, b):
        """
        Define a distance function:
        1) enforce pulled down vertices are labeled with a label which they had before pulling down
        2) cost migrations to be edge independent unit cost
        """
        if a is None or b is None:
            return 0

        if a == b:
            return 0
        else:
            return 1

    def label_f(self, a):
        """
        Determine the label of a given vertex, if it has a unique assignment return it
        else declare it None to be later labeled

        :param self: multi-labeled tree
        :param a: vertex
        """

        label = self.labeling.get(a, [])

        if len(label) == 1:
            return label[0]
        else:
            # this would be the case if it was augmented or unlabeled form the beginning
            return None


def process_model(model, network, character_set, vertex_labeling={}):
    edgelist = []
    if len(vertex_labeling.keys()) == 0:
        # compute (an) optimal vertex label
        vertex_labeling = {}
        for u, v in network.edges:
            if model.reticulations[u, v]() == 1:
                edgelist.append((u, v))
            for c1, c2 in itertools.product(character_set, character_set):
                if (
                    np.abs(model.decisions[u, v, c1, c2]()) > 1e-4
                    and np.abs(model.reticulations[u, v]()) == 1
                ):
                    vertex_labeling[u] = c1
                    vertex_labeling[v] = c2

    # compute the new migration multi-graph
    migration_graph = defaultdict(int)
    for u, v in network.edges:
        if model.reticulations[u, v]() == 1:
            if u not in vertex_labeling or v not in vertex_labeling:
                continue

            if vertex_labeling[u] != vertex_labeling[v]:
                migration_graph[(vertex_labeling[u], vertex_labeling[v])] += 1

    return migration_graph, vertex_labeling, edgelist


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


def resolution_and_labeling_tlp(
    N: nx.DiGraph,
    root,
    alphabet,
    restricted_edges,
    search_complexes: list[tuple[list, list]],
    label_f,
    root_label=None,
):
    reticulation_nodes = set([n for n in N.nodes if N.in_degree(n) > 1])
    reticulation_edges = set([(u, v) for u, v in N.edges if N.in_degree(v) > 1])
    migration_edges = [(i, j) for i in alphabet for j in alphabet if i != j]

    N.add_node("dummy_root")
    N.add_edge("dummy_root", root)

    model = pyo.ConcreteModel()

    # keep things simple in the implementation but in theory you can relax these
    # continuous_edges = [(u,v) for u,v in N.edges if (u,v) not in binary_edges]
    # model.decisions_bin = pyo.Var(binary_edges, alphabet, alphabet, domain=pyo.Binary)

    # x_{u,v,c,c'} for all e=(u,v) in Network, c, c' in alpabet
    model.decisions = pyo.Var(N.edges, alphabet, alphabet, domain=pyo.Binary)
    model.reticulations = pyo.Var(N.edges, domain=pyo.Binary)
    model.migrations = pyo.Var(migration_edges, domain=pyo.Binary, initialize=1)

    # fix all non-reticulation edges
    for u, v in N.edges:
        if (u, v) in reticulation_edges:
            continue
        model.reticulations[u, v].fix(1)

    # maybe fix all non-allowable transitions within binary_edges

    model.edge_constraints = pyo.ConstraintList()
    model.leaf_constraints = pyo.ConstraintList()
    model.interal_vertex_constraints = pyo.ConstraintList()
    model.reticulation_constraints = pyo.ConstraintList()
    model.reticulation_edge_constraints = pyo.ConstraintList()
    model.search_complex_constraints = pyo.ConstraintList()
    model.migration_constraints = pyo.ConstraintList()
    model.root_constraints = pyo.ConstraintList()

    # set the root label to be correct if specified
    if root_label is not None:
        model.root_constraint = pyo.ConstraintList()
        model.root_constraint.add(
            sum(model.decisions["dummy_root", root, c, root_label] for c in alphabet)
            == 1
        )

    # add migration constraints
    for u, v, c1, c2 in model.decisions:
        if c1 == c2:
            continue
        model.migration_constraints.add(
            model.decisions[u, v, c1, c2] <= model.migrations[c1, c2]
        )

    # set reticulation constraints
    for v in reticulation_nodes:
        model.reticulation_constraints.add(
            sum(model.reticulations[(u, v)] for u in N.predecessors(v)) == 1
        )

    # set reticulation and edge constraints
    for u, v in N.edges:
        model.reticulation_edge_constraints.add(
            sum(model.decisions[u, v, i, j] for i in alphabet for j in alphabet)
            == model.reticulations[u, v]
        )

    # set the flow constraints
    for u, v in N.edges:
        for i in alphabet:
            for w in N[v]:
                model.edge_constraints.add(
                    sum(model.decisions[u, v, j, i] for j in alphabet)
                    - sum(model.decisions[v, w, i, j] for j in alphabet)
                    <= 2 - model.reticulations[u, v] - model.reticulations[v, w]
                )

                model.edge_constraints.add(
                    sum(model.decisions[u, v, j, i] for j in alphabet)
                    - sum(model.decisions[v, w, i, j] for j in alphabet)
                    >= -2 + model.reticulations[u, v] + model.reticulations[v, w]
                )

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), if c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), if v is a leaf and c' is not the label of v
    for v in N.nodes:
        # if v is labeled by one thing, we want to disallow any transition ending at v without the correct label
        if not label_f(v):
            continue
        # fix predecessors
        for u in N.predecessors(v):
            for i in alphabet:
                for j in alphabet:
                    if j == label_f(v):
                        continue
                    model.decisions[u, v, i, j].fix(0)
        # fix successors
        for u in N.successors(v):
            for i in alphabet:
                for j in alphabet:
                    if j == label_f(v):
                        continue
                    model.decisions[v, u, j, i].fix(0)

    # search complex constraints
    for nodes, labels in search_complexes:
        complement = [label for label in alphabet if label not in labels]
        # edges between nodes in search complex cannot be the same label
        # cannot be labels outside of the label set
        for i, node in enumerate(nodes):
            for u in nodes[i + 1 :]:
                for label in labels:
                    model.decisions[node, u, label, label].fix(0)
            for u in N.predecessors(node):
                for label in complement:
                    for l in alphabet:
                        model.decisions[u, node, l, label].fix(0)
            for u in N.successors(node):
                for label in complement:
                    for l in alphabet:
                        model.decisions[node, u, label, l].fix(0)

        for label in labels:
            model.search_complex_constraints.add(
                sum(
                    model.decisions[parent, node, j, label]
                    for node in nodes
                    for parent in N.predecessors(node)
                    for j in alphabet
                )
                == 1
            )

    weight_f = lambda i, j: 1 if i != j else 0

    # set objective num_edges + num_migrations
    migration_pattern_parsimony = sum(
        model.migrations[(i, j)] for i in alphabet for j in alphabet if i != j
    )
    migration_graph_parsimony = sum(
        model.decisions[u, v, i, j] * weight_f(i, j)
        for u, v in N.edges
        for i in alphabet
        for j in alphabet
    )

    model.objective_expr = migration_graph_parsimony + migration_pattern_parsimony
    # model.objective_expr = migration_graph_parsimony
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
    return model


def parse_args():
    parser = argparse.ArgumentParser(
        description="Softwired parsimony problem solver using TLP"
    )
    parser.add_argument("network", help="Phylogenetic network in edgelist format")
    parser.add_argument("labeling", help="Leaf sequences in CSV format")
    parser.add_argument("-o", "--output", help="Output prefix")
    parser.add_argument("-r", "--root", help="Root label", default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    try:
        network = nx.read_edgelist(
            args.network, create_using=nx.DiGraph(), data=(("weight", float),)
        )
    except Exception as e:
        network = nx.read_edgelist(args.network, create_using=nx.DiGraph())

    (labeling, alphabet) = parse_multilabelings(args.labeling)
    multi_labeled_tree = MultiLabeledTree.from_tree(network, labeling, alphabet)

    (N, augmented_labeling, restricted_edges, search_complexes) = (
        multi_labeled_tree.augment_tree()
    )

    root = [n for n in N.nodes if N.in_degree(n) == 0]
    if len(root) != 1:
        logger.error("Network must have a single root")
        sys.exit(1)

    root = root[0]

    if not nx.is_directed_acyclic_graph(N):
        logger.error("Network must be a directed acyclic graph.")
        sys.exit(1)

    model = resolution_and_labeling_tlp(
        N,
        root,
        alphabet,
        restricted_edges,
        search_complexes,
        multi_labeled_tree.label_f,
        args.root,
    )

    solver = pyo.SolverFactory("gurobi")
    results = solver.solve(model, tee=True)

    results_dict = {}
    results_dict["objective"] = model.objective()
    results_dict["runtime"] = results.solver.wall_time

    with open(args.output + "_results.json", "w") as f:
        json.dump(results_dict, f)

    (migration_graph, vertex_labeling, edgelist) = process_model(model, N, alphabet)
    write_results(migration_graph, vertex_labeling, args.output)

    with open(f"{args.output}_apr.edgelist", "w+") as f:
        for u, v in edgelist:
            f.write(f"{u}\t{v}\n")
