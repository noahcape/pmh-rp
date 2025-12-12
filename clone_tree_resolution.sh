#!/usr/bin/env bash

TYPE=$1
SUBTYPE=$2

# shift to take additional argument
shift 2

LOC=./data/clone_trees/${TYPE}/${SUBTYPE}/${SUBTYPE}

# solve tlp
gtime -v \
    python \
    ./scripts/tlp.py \
    fast_machina \
    ${LOC}.edgelist \
    ${LOC}.labeling \
    -e 1 -n 1 \
    $@ \
    -o $LOC \
    > ${LOC}_timing.txt 2>&1 

# draw tree
python \
    ./scripts/plots/draw_colored_tree.py \
    ${LOC}_pulled_down.edgelist \
    ${LOC}_vertex_labeling.csv \
    -f "edgelist" \
    -o $LOC

# draw results
dot -Tpng ${LOC}_color_graph.dot > ${LOC}_color_graph.png
dot -Tpng ${LOC}_colored_tree.dot > ${LOC}_colored_tree.png

rm ${LOC}_*.dot


