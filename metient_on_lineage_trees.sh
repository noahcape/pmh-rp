#!/bin/bash
# Run tlp, mach2 and metient on real data

working_dir=/n/fs/ragr-research/projects/pmh-rp
metient=${working_dir}/scripts/processing/run_metient_modified.py
tlp=${working_dir}/scripts/tlp.py
process_metient=${working_dir}/scripts/processing/process_metient_output.py
draw_colored_tree=/n/fs/ragr-research/projects/pmh-rp/scripts/plots/draw_colored_tree.py

data_dir="${working_dir}/laml"

# clones=("68" "28" "30" "17" "20" "21")
clones=("30")
modes=("tlp_normal" "tlp_regularized" "tlp_dag" "tlp_tree")
#todo run mach2
# clones=("28")

echo "Running on $(hostname)"

into_mach2_input_labeling() {
    local i=$1
    local o=$2

    # Convert CSV to TSV and remove the header
    tail -n +2 $i | tr ',' '\t' > $o
}

into_mach2_input_edgelist() {
    local i=$1
    local o=$2
    cut -f1,2 $i > $o
}

# mamba activate metient_gpu
# mamba activate mach2
for clone in ${clones[@]}; do
    # for mode in ${modes[@]}; do
    dir="${data_dir}/m5k_lg${clone}"
    leaf_labeling="${dir}/laml_leaf_labeling.csv"
    edgelist="${dir}/laml_tree_edgelist.tsv"

    echo "Processing: $dir\n"

    metient_timing=${dir}/metient_timing.txt
    metient_out="${dir}/metient"

    tlp_base=$mode
    # mkdir ${dir}/mach2/figures
    out=${dir}/mach2/figures/mach2_0

    # if [[ -f "${dir}/mach2/LL-T-00.location.labeling" ]]; then
    mamba activate tlp
    python \
        $draw_colored_tree \
        $edgelist \
        "${dir}/mach2/LL-T-00.location.labeling" \
        -o $out \
        -f "edgelist" \
        -m
        
    python \
        $draw_colored_tree \
        $edgelist \
        "${dir}/tlp_regularized_exact_vertex_labeling.csv" \
        -o ${dir}/tlp_regularized \
        -f "edgelist" \
        -m
    python \
        $draw_colored_tree \
        $edgelist \
        "${dir}/tlp_normal_vertex_labeling.csv" \
        -o ${dir}/tlp_normal \
        -f "edgelist" \
        -m
    # else
    #     python \
    #         $draw_colored_tree \
    #         $edgelist \
    #         "${dir}/mach2/LL-T-0.location.labeling" \
    #         -o $out \
    #         -f "edgelist" \
    #         -m
    # fi
        

    # dot -Tpng -Gdpi=300 "${out}_color_graph.dot" -o ${out}_color_graph.png
    # dot -Tpng -Gdpi=300 "${out}_colored_tree.dot" -o ${out}_color_tree.png

    # mach2_tree=${dir}/laml.tree
    # mach2_labeling=${dir}/laml.labeling
    # mach2_dir=${dir}/mach2
    # mach2_timing=${dir}/mach2_timing.txt
    # # # mach2
    # into_mach2_input_edgelist $edgelist $mach2_tree 
    # into_mach2_input_labeling $leaf_labeling $mach2_labeling
    # # # solve mach2

    # /usr/bin/time -v \
    #     mach2 \
    #     $mach2_tree \
    #     $mach2_labeling \
    #     -o $mach2_dir \
    #     -p "LL" \
    #     --max_solutions 10 \
    #     > $mach2_timing 2>&1 

    # timeout -k 30s 10h \
    #     /usr/bin/time -v \
    #     python \
    #     $metient \
    #     $edgelist \
    #     $leaf_labeling \
    #     -m $metient_out \
    #     -o $metient_out \
    #     -r "LL" \
    #     > $metient_timing 2>&1 

    # convert the metient tree to a tsv
    # awk '{print "s"$1 "\t" "s"$2}' ${metient_out}_metient_tree.txt > ${metient_out}_metient_tree.tsv

    mamba activate metient_gpu
    python \
        $process_metient \
        $metient_out \
        $metient_timing \
        ${metient_out}_metient_tree.tsv \
        $edgelist \
        "empty"
    # done
done