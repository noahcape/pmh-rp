#!/bin/bash
# Run tlp, mach2 and metient on real data

working_dir=/n/fs/ragr-research/projects/pmh-rp
metient=${working_dir}/scripts/processing/run_metient_modified.py
tlp=${working_dir}/scripts/tlp.py

data_dir="${working_dir}/laml"

for data in "${data_dir}"/*; do
    leaf_labeling="${data}/laml_leaf_labeling.csv"
    edgelist="${data}/laml_tree_edgelist.tsv"

    # skip if LL not found
    if ! grep -q 'LL' "$leaf_labeling"; then
        echo "LL not found in $leaf_labeling — skipping"
        continue
    fi

    echo "Processing: $data - tlp\n"

    # I dont care about timing here but may be useful to compare
    # mach2_timing=${data}/mach2_timing.txt
    # metient_timing=${data}/metient_timing.txt
    # tlp_regularized_timing=${data}/tlp_regularized_timing.txt
    # tlp_normal_timing=${data}/tlp_normal_timing.txt
    # tlp_dag_timing=${data}/tlp_dag_timing.txt
    # tlp_tree_timing=${data}/tlp_tree_timing.txt

    #out files
    metient_out="${data}/metient"
    tlp_regularized_out=${data}/tlp_regularized
    tlp_normal_out=${data}/tlp_normal
    tlp_dag_out=${data}/tlp_dag
    tlp_tree_out=${data}/tlp_tree


    # mach2
    # into_mach2_input_edgelist $edgelist $mach2_tree 
    # into_mach2_input_labeling $leaf_labeling $mach2_labeling

    # solve mach2
    # /usr/bin/time -v \
    #     conda run -n mach2 \
    #     mach2 \
    #     $mach2_tree \
    #     $mach2_labeling \
    #     -o $mach2_dir \
    #     -p "LL" \
    #     --max_solutions 10 \
    #     > $mach2_timing 2>&1 

    # metient with 3h timeout
    # timeout -k 30s 3h \
    #     /usr/bin/time -v \
    #     env CUDA_VISIBLE_DEVICES="" \
    #     conda run -n metient \
    #     python \
    #     $metient \
    #     $edgelist \
    #     $leaf_labeling \
    #     -o $metient_out \
    #     -r "LL" \
    #     > $metient_timing 2>&1 

    # solve tlp - normal
    # /usr/bin/time -v \
    #     conda run -n tlp \
    #     python \
    #     $tlp \
    #     fast_machina \
    #     $edgelist \
    #     $leaf_labeling \
    #     -o $tlp_normal_out \
    #     -l "LL" \
    #     > $tlp_normal_timing 2>&1 

    # solve tlp - edge regularized normal
    /usr/bin/time -v \
        conda run -n tlp \
        python \
        $tlp \
        fast_machina \
        $edgelist \
        $leaf_labeling \
        -o $tlp_regularized_out \
        -e 1 \
        -l "LL" \
        > $tlp_regularized_timing 2>&1 

    # # solve tlp - tree constrained
    # /usr/bin/time -v \
    #     conda run -n tlp \
    #     python \
    #     $tlp \
    #     fast_machina \
    #     $edgelist \
    #     $leaf_labeling \
    #     -o $tlp_tree_out \
    #     -c "polyclonal_tree" \
    #     -l "LL" \
    #     > $tlp_tree_timing 2>&1 

    # # solve tlp - dag constrained
    # /usr/bin/time -v \
    #     conda run -n tlp \
    #     python \
    #     $tlp \
    #     fast_machina \
    #     $edgelist \
    #     $leaf_labeling \
    #     -o $tlp_dag_out \
    #     -c "polyclonal_dag" \
    #     -l "LL" \
    #     > $tlp_dag_timing 2>&1 
done