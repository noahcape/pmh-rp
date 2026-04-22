#!/bin/bash

working_dir=/n/fs/ragr-research/projects/pmh-rp/clone_trees
tlp=/n/fs/ragr-research/projects/pmh-rp/scripts/tlp.py
clone_tree_tlp=/n/fs/ragr-research/projects/pmh-rp/scripts/dag_resolution_labeling.py
metient=/n/fs/ragr-research/projects/pmh-rp/scripts/processing/run_metient_clone_trees.py
process_metient=/n/fs/ragr-research/projects/pmh-rp/scripts/processing/process_metient_output.py
draw_colored_tree=/n/fs/ragr-research/projects/pmh-rp/scripts/plots/draw_colored_tree.py
data=("A10" "A22")


for dir in ${data[@]}; do
    data_dir="${working_dir}/${dir}"
    tree="${data_dir}/${dir}.tree"
    labeling="${data_dir}/${dir}.observed.labeling"
    mach2_dir="${data_dir}/mach2"
    metient_dir="${data_dir}/metient"
    timing_file="${data_dir}/mach2_timing.txt"
    tlp_timing_file="${data_dir}/tlp_clone_tree_timing.txt"
    metient_timing_file="${data_dir}/metient_timing.txt"


    # find root

    # tlp
    mamba activate tlp
    /usr/bin/time -v \
        python \
        $clone_tree_tlp \
        $tree \
        $labeling \
        -o $tlp_out \
        -r $root \
        > $tlp_timing_file 2>&1


    mamba activate metient_gpu
    /usr/bin/time -v \
        python \
        $metient \
        $tree \
        $labeling \
        -x 1 \
        -r $root \
        -m $metient_dir \
        -o $data_dir \
        > $metient_timing_file 2>&1

    # convert the metient tree to a tsv
    awk '{print "s"$1 "\t" "s"$2}' ${data_dir}/metient_tree.txt > ${data_dir}/metient_tree.tsv

    # metient
    python \
        $process_metient \
        $metient_dir \
        $metient_timing_file \
        ${data_dir}/metient_tree.tsv \
        $tree \
        "empty" \
        >> $metient_timing_file


    # solve mach2
    mamba deactivate tlp
    mamba activate mach2
    /usr/bin/time -v \
        conda run -n mach2 \
        mach2 \
        $tree \
        $labeling \
        -o $mach2_dir \
        -p $root \
        --max_solutions 10 \
        > $timing_file 2>&1 


done
