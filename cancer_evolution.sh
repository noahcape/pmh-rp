#!/usr/bin/env bash

SUBTYPE=$1

# shift to take additional argument
shift 1

LOC=./examples/cancer_evolution/${SUBTYPE}/${SUBTYPE}

# solve tlp
gtime -v \
    python \
    ./scripts/tlp.py \
    fast_machina \
    ${LOC}_tree_edgelist.tsv \
    ${LOC}_leaf_labeling.csv \
    -e 1 \
    $@ \
    -o $LOC \
    > ${LOC}_timing.txt 2>&1 

# draw tree
python \
    ./scripts/plots/draw_colored_tree.py \
    ${LOC}_tree_edgelist.tsv \
    ${LOC}_vertex_labeling.csv \
    -f "edgelist" \
    -o $LOC

# draw results
dot -Tpng ${LOC}_color_graph.dot > ${LOC}_color_graph.png
dot -Tpng ${LOC}_colored_tree.dot > ${LOC}_colored_tree.png

rm ${LOC}_*.dot


