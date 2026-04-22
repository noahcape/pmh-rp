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

from dag_resolution_labeling import MultiLabeledTree
from dag_resolution_labeling import resolution_and_labeling_tlp

def test_big():
    edgelist = [
        ("grey", "red"),
        ("grey", "dark_brownC"),
        ("red", "mid_blue"),
        ("red", "light_green"),
        ("light_green", "dark_purple"),
        ("light_green", "dark_blue"),
        ("dark_purple", "orange"),
        ("dark_purple", "mid_green"),
        ("dark_purple", "dark_green"),
        ("orange", "pink"),
        ("orange", "light_purpleE"),
        ("mid_green", "goldF"),
        ("dark_green", "light_blue"),
        ("dark_green", "yellow"),
        ("dark_green", "light_brownI"),
    ]

    labeling = {
        "grey": ["seminal_vesicle", "prostate"],
        "red": ["left_pelvic_lymph_node_5", "seminal_vesicle", "prostate"],
        "mid_blue": ["bladder", "right_pelvic_lymph_node_12", "seminal_vesicle"],
        "light_green": [
            "pelvic_lymph_node_7",
            "right_pelvic_lymph_node_12",
            "left_pelvic_lymph_node_8",
            "left_pelvic_lymph_node_5",
            "right_adrenal_gland",
            "bladder",
        ],
        "dark_purple": [
            "pelvic_lymph_node_7",
            "right_pelvic_lymph_node_12",
            "left_pelvic_lymph_node_8",
            "left_adrenal_gland",
            "left_pelvic_lymph_node_5",
            "right_adrenal_gland",
            "bladder",
        ],
        "dark_blue": ["bladder", "right_pelvic_lymph_node_12", "seminal_vesicle"],
        "orange": [
            "left_humerus_bone_marrow",
            "right_pelvic_lymph_node_12",
            "bladder",
            "left_adrenal_gland",
        ],
        "mid_green": ["right_adrenal_gland", "left_adrenal_gland"],
        "dark_green": ["left_pelvic_lymph_node_5", "left_pelvic_lymph_node_8"],
        "light_blue": ["left_pelvic_lymph_node_5", "pelvic_lymph_node_7"],
        "yellow": ["left_pelvic_lymph_node_5", "pelvic_lymph_node_7"],
        "pink": ["left_humerus_bone_marrow", "right_pelvic_lymph_node_12"],
        "dark_brownC": ["prostate"],
        "goldF": ["right_adrenal_gland"],
        "light_brownI": ["left_pelvic_lymph_node_8"],
        "light_purpleE": ["left_adrenal_gland"],
    }

    character_set = list(set([item for value in labeling.values() for item in value]))

    clone_tree = MultiLabeledTree.from_edgelist(
        edgelist, labeling, character_set
    )

    (N, augmented_labeling, restricted_edges, search_complexes) = clone_tree.augment_tree()
    # for u,v in N.edges:
    #     print(u,v)
    # for u in N.nodes:
    #     if u in augmented_labeling.keys():
    #         print(u, augmented_labeling[u])

    root = [n for n in N.nodes if N.in_degree(n) == 0]
    if len(root) != 1:
        print("Network must have a single root")
        exit(1)

    root = root[0]

    if not nx.is_directed_acyclic_graph(N):
        logger.error("Network must be a directed acyclic graph.")
        sys.exit(1)

    alphabet = clone_tree.character_set
    model = resolution_and_labeling_tlp(N, root, alphabet, restricted_edges, search_complexes, clone_tree.label_f)

    solver = pyo.SolverFactory('gurobi')
    results = solver.solve(model, tee=True)

    for u, v in N.edges:
        for i in alphabet:
            for j in alphabet:
                if np.abs(model.decisions[u, v, i, j]()) > 1e-4:
                    print(u, "labeled: ", i, v, "labeled: ", j)

    results_dict = {}
    results_dict['objective'] = model.objective()
    results_dict['runtime'] = results.solver.wall_time
    print(results_dict)

    
if __name__ == "__main__":

    edgelist = [
        (1,2),
        (2,5),
        (5,3),
        (5,4),
        # (3,5),
        # (5,6),
        # (6,8)
        # (4,7),
    ]

    labeling = {
        # 7: ["A"],
        1: ["A"],
        2: ["B"],
        5: ["B", "C"],
        # 2: ["A", "B", "C", "D", "E"],
        # 3: ["A", "B"],
        3: ["A", "B"],
        4: ["B", "C"],
    }

    character_set = ["A", "B", "C"]
    clone_tree = MultiLabeledTree.from_edgelist(
        edgelist, labeling, ["A", "B", "C"]
    )

    (N, augmented_labeling, restricted_edges, search_complexes) = clone_tree.augment_tree()

    for u,v in N.edges:
        print(u,v)
    for u in N.nodes:
        if u in augmented_labeling.keys():
            print(u, augmented_labeling[u])

    root = [n for n in N.nodes if N.in_degree(n) == 0]
    if len(root) != 1:
        print("Network must have a single root")
        exit(1)

    root = root[0]


    if not nx.is_directed_acyclic_graph(N):
        logger.error("Network must be a directed acyclic graph.")
        sys.exit(1)

    alphabet = clone_tree.character_set
    model = resolution_and_labeling_tlp(N, root, alphabet, restricted_edges, search_complexes, clone_tree.label_f)

    solver = pyo.SolverFactory('gurobi')
    results = solver.solve(model, tee=True)

# def process_model(model, tree, character_set, vertex_labeling={}):
    vertex_labeling = {}
    if len(vertex_labeling.keys()) == 0:
        # compute (an) optimal vertex label
        vertex_labeling = {}
        for u,v in N.edges:
            if model.reticulations[u,v]() == 1:
                print(u, v)
            for c1, c2 in itertools.product(character_set, character_set):
                if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                    vertex_labeling[u] = c1
                    vertex_labeling[v] = c2

    # compute the new migration multi-graph
    migration_graph = defaultdict(int)
    for u, v in N.edges:
        if u not in vertex_labeling or v not in vertex_labeling:
            continue

        if vertex_labeling[u] != vertex_labeling[v]:
            migration_graph[(vertex_labeling[u], vertex_labeling[v])] += 1

    # print(migration_graph)
    print(vertex_labeling)

    results_dict = {}
    results_dict['objective'] = model.objective()
    results_dict['runtime'] = results.solver.wall_time
    print(results_dict)
