import sys
import argparse
import itertools
import json

import pandas as pd
import numpy as np
import networkx as nx
import pyomo.environ as pyo
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict
from gurobipy import GRB

import matplotlib.pyplot as plt
from kneed import KneeLocator

def find_knee(x, y, s, curve):
    """
    Locate the knee of a concave increasing function to identify the optimal parameters
    
    :param x: regularizing parameter (l)
    :param y: optimal value = migrations + l * edges
    :param s: sensitivity parameter for knee location

    Return the knee (x, y)
    """

    # TODO: should be concave but make sure - plot this
    # how do we define a pareto front - just the sum of the two values or consider the two values seperately
    knee = KneeLocator(x, y, s, curve=curve, direction="increasing")

    return (knee.knee, knee.knee_y)

def process_model(model, tree, character_set, vertex_labeling={}):
    if len(vertex_labeling.keys()) == 0:
        # compute (an) optimal vertex label
        vertex_labeling = {}
        for u,v in tree.edges:
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

    print(f"<------ Migrations in labeling: {migrations} ------>\n<------ Number of edges in migration graph: {migation_pattern_edges} ------>")


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
                raise ValueError(f"Lines in labeling file for clone tree ancestral reconstruction expected to be of the form <node> <label> (<label>)*")

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

    def pull_down(self, out=None):
        """
        Pull down all multilabeled nodes. Create a dictionary associating each
        pulled down node with the set of labeled with which is can be labeled with.
        """

        # list of nodes along with the set of labeled which they should be constrained to be labeled by
        # self.pulled_down_vertices = {}

        new_edges = []
        
        
        for v in self.tree.nodes:
            labels = self.labeling.get(v)

            if labels is not None and len(labels) > 1:
                # for each label create a new node
                for l in labels:
                    new_node = f"{v}_{l}"
                    # add new labeled vertex to the labeling
                    self.labeling[new_node] = [l]
                    
                    # keep track of new nodes to add
                    new_edges.append((v,new_node))

        self.tree.add_edges_from(new_edges)

        if out is not None:
            # here we need to print out the new edge list
            with open(f"{out}_pulled_down.edgelist", "w+") as f:
                for (u,v) in self.tree.edges:
                    f.write(f"{u}\t{v}\n")


    def dist_f(self, t, a, b):
        """
        Define a distance function:
        1) enforce pulled down vertices are labeled with a label which they had before pulling down
        2) cost migrations to be edge independent unit cost
        """
        if a is None or b is None:
            return 0
        
        (u,v) = t
        valid_u_labels = self.labeling.get(u, [])
        valid_v_labels = self.labeling.get(v, [])

        # no penalty as long as the node has no restrictions or if the label is nodes restrictions
        def no_penalty(restricted_list, label):
            return len(restricted_list) == 0 or label in restricted_list
        
        # if either node requires a penalty then dist approx infinity
        if no_penalty(valid_u_labels, a) and no_penalty(valid_v_labels, b):
            if a == b:
                return 0
            else:
                return 1
        
        # TODO: set to inf then set x_uvab = 0      
        return len(self.tree.nodes)
 
        
    def label_f(self, a):
        label = self.labeling.get(a, [])

        if len(label) == 0 or len(label) > 1:
            return None
        else:
            return label[0]



    def print_tree(self):
        print(self.tree.edges)
        print(self.labeling)

    # TODO: just use normal pareto_front_solver
    def solve_tlp(self, l):
        solver = pyo.SolverFactory('gurobi_persistent')

        root = [v for v in self.tree.nodes if self.tree.in_degree(v) == 0][0]
        
        model = create_tree_labeling_polytope(self.tree, root, self.character_set, self.leaf_f, self.dist_f)
        model = set_regularized_parsimony_objective(self.tree, self.character_set, model, self.dist_f, l)

        solver.set_instance(model)
        logger.info("Solving the full MILP model")
        solver.options["MIPGap"] = 0.15
        solver.solve(model, tee=True, warmstart=True)

        return model

def solve_pareto_front(tree, root, character_set, dist_f, leaf_f, root_label, output, weights, save_fig=True, mip_gap=0):
    """
    Solve the tlp for different weights on the number of edges for the regularized objective function
    Warm start with the solution to the previous instance as with increasing lambda it is a lower bound on the solution

    
    :param tree: tree topology of problem instance
    :param root: root of the tree
    :param character_set: set of labels
    :param dist_f: distance function for scoring parsimony
    :param leaf_f: function to get label of a node if it is a leaf
    :param root_label: possible label of the root node
    :param output: output path
    :param weights: weights to use for weighting edge number regularization
    """

    y_parsimony_score = []
    y_migration_pattern_edges = []
    x_regularization_weights = weights
    solutions = {}

    # solve first with the first weight
    solver = pyo.SolverFactory('gurobi_persistent')
    model = create_tree_labeling_polytope(tree, root, character_set, leaf_f, dist_f, root_label)
    model = set_regularized_parsimony_objective(tree, character_set, model, dist_f, weights[0])
    
    solver.set_instance(model)
    logger.info("Solving the full MILP model")
    solver.options["MIPGap"] = mip_gap
    solver.solve(model, tee=True, warmstart=True)

    (migration_graph, vertex_labeling) = process_model(model, tree, character_set)
    migrations = sum(migration_graph.values())
    solutions[weights[0]] = (migration_graph, vertex_labeling) 
    y_parsimony_score.append(migrations)
    y_migration_pattern_edges.append(len(migration_graph.keys()))

    # then again with new weights with warm start of previous solution
    for lambda_ in weights[1:]:
        # delete previous obj
        model.del_component(model.objective_expr)
        model.del_component(model.objective)

        # TODO: ensure that we dont have the create a new solver instance
        # solver = pyo.SolverFactory('gurobi_persistent')
        
        # set a new objective
        model = set_regularized_parsimony_objective(tree, character_set, model, dist_f, lambda_, False)
        solver.set_instance(model)

        logger.info("Solving the full MILP model")
        solver.options["MIPGap"] = mip_gap
        solver.solve(model, tee=True, warmstart=True)

        (migration_graph, vertex_labeling) = process_model(model, tree, character_set)
        migrations = sum(migration_graph.values())
        solutions[lambda_] = (migration_graph, vertex_labeling) 
        y_parsimony_score.append(migrations)
        y_migration_pattern_edges.append(len(migration_graph.keys()))

    print(x_regularization_weights)
    print(y_parsimony_score)
    print(y_migration_pattern_edges)
    # find the best solution on the pareto front
    try:
        (knee_x, knee_y) = find_knee(x_regularization_weights, y_parsimony_score, 0.5, "concave")
        if knee_x == None or knee_y == None:
            (knee_x, knee_y) = find_knee(x_regularization_weights, y_parsimony_score,  0.5, "convex")
    except Exception as _:
        knee_y = 1

    if knee_x == None or knee_y == None:
        knee_x = 1
    
    # write out the pareto front
    if save_fig:
        plt.plot(x_regularization_weights, y_parsimony_score, marker='o', linewidth=2, markersize=1)
        plt.ylabel("Migrations")
        plt.plot([knee_x], [knee_y], marker='*', color ='red', markersize=10)
        plt.xlabel("Regularization Weight")
        plt.title(f"Regularized Parsimony Pareto Front (s=0.5)")
        plt.savefig(f"{output}_pareto_front.png")
        plt.clf()
    
    (migration_graph, vertex_labeling) = solutions[knee_x]
    write_results(migration_graph, vertex_labeling, output, { "regularizer_weight": knee_x })


def is_leaf(T, node):
    return len(T[node]) == 0


def set_parsimony_objective(T, character_set, model, dist_f):
    """
    Set the objective of the model to be the number of migrations in the tree.
    """

    # \sum_{u,v} \sum_{c,c'} x_{u,v,c,c'} * d((u,v), c,c')
    migrations = sum(model.decisions[u,v,c,c2] * dist_f((u,v), c, c2) for u,v in T.edges for c in character_set for c2 in character_set)

    model.objective_expr = migrations
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
    return model

def set_regularized_parsimony_objective(T, character_set, model, dist_f, regularization_param, add_migrations=True):
    """
    Set the objective of the model to be the sum of the number of edges in the migration pattern
    and the number of migrations in the tree.
    """
    logger.info(f"Setting object with weight: {regularization_param}")

    if add_migrations:
        # append migrations since we need it in our objective
        append_migrations(model, tree, character_set)

    # set objective num_edges + num_migrations
    num_edges = sum(model.migrations[(i,j)] for i in character_set for j in character_set if i != j)

    # \sum_{u,v} \sum_{c,c'} x_{u,v,c,c'} * d((u,v), c,c')
    migrations = sum(model.decisions[u,v,c,c2] * dist_f((u,v), c, c2) for u,v in T.edges for c in character_set for c2 in character_set)

    model.objective_expr = (regularization_param * num_edges) + migrations
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
    return model

"""
Creates the tree labeling polytope whose vertices correspond
to the set of solutions to the (unconstrained) maximum parsimony 
problem.
"""
def create_tree_labeling_polytope(T, root, character_set, label_f, dist_f, root_label=None):
    model = pyo.ConcreteModel()

    logger.info(f"Creating tree labeling polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")
    logger.info(f"Character set has size: {len(character_set)}")

    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()
    # set domain to be [0, 1]

    # model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.Binary, initialize=1)
    model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0)

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w,c
    model.edge_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                model.edge_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set) - sum(model.decisions[v, w, c, c2] for c2 in character_set) == 0
                )

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), v is a leaf, and c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), v is a leaf, and c' is not the label of v
    model.leaf_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        if label_f(v) == None: continue
        for c in character_set:
            if c == label_f(v):
                model.leaf_constraints.add(sum(model.decisions[u, v, c2, c] for c2 in character_set) == 1)
            else:
                model.leaf_constraints.add(sum(model.decisions[u, v, c2, c] for c2 in character_set) == 0)

    # set root label if provided
    if root_label is not None:
        model.root_constraint = pyo.ConstraintList()
        model.root_constraint.add(sum(model.decisions["dummy_root", root, c, root_label] for c in character_set) == 1)

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
        if c1 == c2: continue
        model.migration_constraints.add(
            model.decisions[u, v, c1, c2] <= model.migrations[c1, c2]
        )

"""
Constrain the tree labeling polytope to only allow certain migration graphs
by adding constraints to the model using a Gurobi callback.
"""
def setup_fast_machina_constraints(solver, model, character_set, constraint_type):
    def dag_callback(model, gb_model, where):
        if where != GRB.Callback.MIPSOL:
            return

        # load solutions
        gb_model.cbGetSolution([model.migrations[i, j] for i in character_set for j in character_set if i != j])

        G = nx.DiGraph()
        for i in character_set:
            G.add_node(i)

        for i in character_set:
            for j in character_set:
                if i == j: continue
                if model.migrations[i, j].value > 0.5:
                    G.add_edge(i, j)

        try:
            S = nx.find_cycle(G, orientation="original")
            S = [i for (i,_,_) in S]

            logger.info(f"Adding constraint for cycle {S}.")

            cons = model.migration_graph_constraints.add(
                sum(model.migrations[(S[i], S[(i+1) % len(S)])] for i in range(len(S))) <= len(S) - 1
            )

            gb_model.cbLazy(cons)
        except nx.exception.NetworkXNoCycle:
            pass

    model.migration_graph_constraints = pyo.ConstraintList()
    if constraint_type == "polyclonal_tree" or constraint_type == "monoclonal_tree":
        S = character_set
        model.migration_graph_constraints.add(
            sum(model.migrations[(i, j)] for i in S for j in S if i != j) == len(S) - 1
        )
        solver.set_instance(model)
    elif constraint_type == "polyclonal_dag" or constraint_type == "monoclonal_dag":
        solver.set_instance(model)
        solver.set_gurobi_param("LazyConstraints", 1)
        solver.set_callback(dag_callback)
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

    return

"""
Solves a generalization of the convex recoloring
problem using the formulation of Campelo et al. 2016.
"""
def campelo_et_al(tree, character_set, leaf_f, dist_f, root, args, integral=True):
    rooted_T = tree.copy() 
    T = rooted_T.to_undirected()

    model = pyo.ConcreteModel()
    if integral:
        model.x = pyo.Var(T.nodes, character_set, domain=pyo.Binary)
    else:
        model.x = pyo.Var(T.nodes, character_set, domain=pyo.NonNegativeReals, bounds=(0, 1))

    logger.info(f"Creating Campelo et al. 2016 formulation for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")

    model.path_constraints = pyo.ConstraintList()
    for u in tqdm(T.nodes):
        for v in T.nodes:
            if u == v: continue
            path = nx.shortest_path(T, u, v)
            for w in path[1:-1]:
                for c in character_set:
                    model.path_constraints.add(model.x[u, c] + model.x[v, c] - model.x[w, c] <= 1)

    for u in T.nodes:
        model.path_constraints.add(
            sum(model.x[u, c] for c in character_set) == 1
        )

    def weight_f(u, c):
        if rooted_T.out_degree(u) != 0: return 0
        if c == leaf_f(u): return 0 
        return 1

    model.objective_expr = sum(weight_f(u, c) * model.x[u, c] for u in T.nodes for c in character_set)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)

    solver = pyo.SolverFactory('gurobi_persistent')
    solver.set_instance(model)
    solver.solve(model, tee=True)

    vertex_labeling = {}
    for u in T.nodes:
        for c in character_set:
            if np.abs(model.x[u, c]()) > 1e-4:
                vertex_labeling[u] = c

    return vertex_labeling, model.objective()

def parsimonious_relabeling(tree, character_set, leaf_f, dist_f, root, args, integral=True):
    T = tree.copy()

    if args.k is None:
        k = len(character_set) - 1

    model = pyo.ConcreteModel()

    logger.info(f"Creating tree labeling polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")
    logger.info(f"Character set has size: {len(character_set)}")

    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    pendant_edges = [(u, v) for u, v in T.edges if is_leaf(T, v)]

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()

    # set domain to be [0, 1]
    if integral:
        model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.Binary, initialize=0)
        model.relabelings = pyo.Var(pendant_edges, domain=pyo.Binary, initialize=0)
    else:
        model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0)
        model.relabelings = pyo.Var(pendant_edges, domain=pyo.NonNegativeReals, bounds=(0, 1))

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w,c
    model.edge_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                model.edge_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set) - sum(model.decisions[v, w, c, c2] for c2 in character_set) == 0
                )

    # require \sum_{c,c'} x_{u,v,c,c'} = 1 for all e=(u,v)
    model.sum_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        model.sum_constraints.add(sum(model.decisions[u, v, c1, c2] for c2 in character_set for c1 in character_set) == 1)

    # require leaves that are not relabeled to have the correct label
    model.leaf_constraints = pyo.ConstraintList()
    for u, v in pendant_edges:
        for c in character_set:
            if c == leaf_f(v):
                model.leaf_constraints.add(sum(model.decisions[u, v, c2, c] for c2 in character_set) == model.relabelings[u, v])

    # c^T x_{u,v,i,j} \leq k
    model.sum_constraints.add(
        sum(model.decisions[u, v, c1, c2] * dist_f((u, v), c1, c2) for u, v in T.edges for c1 in character_set for c2 in character_set) <= k
    )

    for color in character_set:
        leaves = [v for v in T.nodes if is_leaf(T, v) and leaf_f(v) == color]
        if len(leaves) == 0: continue
        model.sum_constraints.add(
            sum(model.relabelings[u, v] for u, v in pendant_edges if v in leaves) >= 1
        )

    model.objective_expr = sum(1 - model.relabelings[u, v] for u, v in pendant_edges)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)

    solver = pyo.SolverFactory('gurobi_persistent')
    solver.set_instance(model)
    solver.solve(model, tee=True, warmstart=True)

    vertex_labeling = {}
    for u, v in T.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2
    
    return vertex_labeling, model.objective()

"""
Solves a generalization of the MACHINA parsimonious migration 
history problem using the tree labeling polytope.
TODO: next step if to simplify the main method just in here add regularization parameter
"""
def fast_machina(tree, character_set, leaf_f, dist_f, root, args, mip_gap=0.15):
    tree = tree.copy()

    """ Step 1: Setup the MILP using Gurobi callbacks, if necessary """
    solver = pyo.SolverFactory('gurobi_persistent')

    if "dummy_root" in tree.nodes:
        tree.remove_node("dummy_root")

    model = create_tree_labeling_polytope(tree, root, character_set, leaf_f, dist_f, root_label=args.label)
    # set objective for min migrations
    model = set_parsimony_objective(tree, character_set, model, dist_f)
    if args.constraints != "none":
        append_migrations(model, tree, character_set)
    if args.constraints.startswith("monoclonal"):
        for c1, c2 in model.migrations:
            model.migration_constraints.add(
                sum(model.decisions[u, v, c1, c2] for u, v in tree.edges) <= model.migrations[c1, c2]
            )

    added_constraints = []
    if args.constraints != "none":
        setup_fast_machina_constraints(solver, model, character_set, args.constraints)
    else:
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

    subparsers = parser.add_subparsers(dest="method", help="Methods")
    subparsers.required = True

    # fastMACHINA subparser
    fast_machina_parser = subparsers.add_parser("fast_machina", help="fastMACHINA")
    fast_machina_parser.add_argument("tree", help="Tree in edgelist format")
    fast_machina_parser.add_argument("labels", help="Leaf labeling as a CSV file")
    fast_machina_parser.add_argument("-c", "--constraints", help="Migration graph constraints",
                                choices=["polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag", "none"],
                                default="none")
    fast_machina_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    fast_machina_parser.add_argument("-r", "--root", help="Root of the tree", default=None)
    fast_machina_parser.add_argument("-l", "--label", help="Root label", default=None)
    fast_machina_parser.add_argument("-e", "--edge-regularized", help="Regularize normal parsimony with number of edges in migraton graph", type=int, default=0)
    fast_machina_parser.add_argument("-w", "--regularizer-weight", help="How much to weigh edge regularization", type=float, default=1.0)
    fast_machina_parser.add_argument("-n", "--clone-tree-instance", help="Expect a clone tree with internal and multilabels", type=int, default=0)

    # parsimonious relabeling subparser
    parsimonious_relabeling = subparsers.add_parser("convex_recoloring", help="parsimoniousRelabeling")
    parsimonious_relabeling.add_argument("tree", help="Tree in edgelist format")
    parsimonious_relabeling.add_argument("labels", help="Leaf labeling as a CSV file")
    parsimonious_relabeling.add_argument("-o", "--output", help="Output prefix", default="result")
    parsimonious_relabeling.add_argument("-r", "--root", help="Root of the tree", default="root")
    parsimonious_relabeling.add_argument("-k", help="Weighted parsimony constraint", default=None, type=float)
    parsimonious_relabeling.add_argument("-w", "--weights", help="Weight of transitioning between labels", default=None)
    parsimonious_relabeling.add_argument("-m", "--mode", help="Mode", choices=["campelo", "tlp"], default="tlp")

    args = parser.parse_args()
    print(args)
    return args

if __name__ == "__main__":
    args = parse_arguments()

    try:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))
    except Exception as e:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph())

    # store these booleans for convenience
    regularize = args.edge_regularized == 1
    clone_tree_instance = args.clone_tree_instance == 1

    if clone_tree_instance:
        # parse labeling and character set a special way for clone tree instances
        (labeling, character_set) = parse_multilabelings(args.labels)
    else:
        # fall back to normal parsing for phylogenetic (non-clone tree) tree instance
        labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")

        # load the character set as the set of all leaf labels
        character_set = labels_csv[labels_csv["label"] != "None"]["label"].unique()

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

    if clone_tree_instance:
        clone_tree = MultiLabeledTree.from_tree(tree, labeling, character_set)
        # pull down tree
        clone_tree.pull_down(args.output)

        # check that all leaves are labeled
        if not all([clone_tree.label_f(node) is not None for node in clone_tree.tree.nodes if is_leaf(clone_tree.tree, node)]):
            unlabeled_leaves = [node for node in clone_tree.tree.nodes if is_leaf(clone_tree.tree, node) and clone_tree.label_f(node) is None]
            raise ValueError(f"Leaves {unlabeled_leaves} are unlabeled.")

        # solve tlp with regularization and get solved model
        # TODO: what weights to use
        weights = [0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1, 1.5, 2]
        root = [v for v in clone_tree.tree.nodes if clone_tree.tree.in_degree(v) == 0][0]

        mip_gap = 0.05 if len(clone_tree.tree.nodes) < 25 else 0.15

        # solve for the pareto front with respect to the regularization parameter
        solve_pareto_front(clone_tree.tree, root, character_set, clone_tree.dist_f, clone_tree.label_f, args.label, args.output, weights, mip_gap=mip_gap)
        
    else:

        # defines the distance function between characters x and y along an edge e
        def dist_f(e, x, y):
            if x is None or y is None:
                return 0

            if x == y:
                return 0

            return 1
        
        def is_int(s):
            try:
                int(s)
                return True
            except ValueError:
                return False

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

        if not all([leaf_f(node) is not None for node in tree.nodes if is_leaf(tree, node)]):
            unlabeled_leaves = [node for node in tree.nodes if is_leaf(tree, node) and leaf_f(node) is None]
            raise ValueError(f"Leaves {unlabeled_leaves} are unlabeled.")

        if args.edge_regularized == 0:
            if len(character_set) == 1:
                logger.warning("Character set has size 1, inferring trivial labeling.")
                vertex_labeling = {node: character_set[0] for node in tree.nodes}
                lp_obj = 0
                obj = 0
            else:
                # computes the vertex labeling using the specified method
                if args.method == "fast_machina":
                    lp_obj = None
                    vertex_labeling, obj = fast_machina(tree, character_set, leaf_f, dist_f, root, args)
                elif args.method == "convex_recoloring":
                    if args.mode == "tlp":
                        _, lp_obj = parsimonious_relabeling(tree, character_set, leaf_f, dist_f, root, args, integral=False)
                        vertex_labeling, obj = parsimonious_relabeling(tree, character_set, leaf_f, dist_f, root, args)
                    else:
                        _, lp_obj = campelo_et_al(tree, character_set, leaf_f, dist_f, root, args, integral=False)
                        vertex_labeling, obj = campelo_et_al(tree, character_set, leaf_f, dist_f, root, args)

            (migration_graph, vertex_labeling) = process_model(None, tree, character_set, vertex_labeling)
            write_results(migration_graph, vertex_labeling, args.output)
        else:
            # solve for the pareto front with respect to the regularization parameter
            # TODO: what weights to use
            weights = [0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1]
            mip_gap = 0.05 if len(tree.nodes) < 25 else 0.15
            solve_pareto_front(tree, root, character_set, dist_f, leaf_f, args.label, args.output, weights, mip_gap=mip_gap)