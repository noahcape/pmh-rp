#!/usr/bin/env bash
# Run tlp, mach2 and metient on real data

working_dir=/n/fs/ragr-research/projects/pmh-rp
metient=${working_dir}/scripts/processing/run_metient_modified.py
tlp=${working_dir}/scripts/tlp.py
process_metient=${working_dir}/scripts/processing/process_metient_output.py

data_dir="${working_dir}/laml"

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

dirs=("m5k_lg68", "m5k_lg17", "m5k_lg30", "m5k_lg28", "m5k_lg20", "m5k_lg21")

for data in "${data_dir}"/*; do
    leaf_labeling="${data}/laml_leaf_labeling.csv"
    edgelist="${data}/laml_tree_edgelist.tsv"

    if [[ "$(basename "$data")" != "m5k_lg68" ]]; then
        continue
    fi

    # metient - takes too long if more than 100
    # skip if file has more than 100 lines
    # if [ "$(wc -l < "$leaf_labeling")" -gt 400 ]; then
    #     continue
    # fi

    # skip if LL not found
    if ! grep -q 'LL' "$leaf_labeling"; then
        # echo "LL not found in $leaf_labeling — skipping"
        continue
    fi


    #timing files
    mach2_timing=${data}/mach2_timing.txt
    metient_timing=${data}/metient_timing.txt
    tlp_regularized_timing=${data}/tlp_regularized_timing_exact.txt
    tlp_normal_timing=${data}/tlp_normal_timing.txt
    tlp_dag_timing=${data}/tlp_dag_timing.txt
    tlp_tree_timing=${data}/tlp_tree_timing.txt

    #out files
    metient_out="${data}/metient"
    tlp_regularized_out=${data}/tlp_regularized_exact
    tlp_normal_out=${data}/tlp_normal
    tlp_dag_out=${data}/tlp_dag
    tlp_tree_out=${data}/tlp_tree


    # mach2_dir=${data}/mach2
    # if [ ! -d $mach2_dir ]; then
    #     mkdir $mach2_dir
    # fi

    # if [ $(ls -A "$mach2_dir" | wc -l) -eq 0 ]; then
    #     echo "Directory is empty ${mach2_dir}"

    #     echo "Processing: $data - mach2\n"

    #     mach2_tree=${data}/laml.tree
    #     mach2_labeling=${data}/laml.labeling

    #     # mach2
    #     into_mach2_input_edgelist $edgelist $mach2_tree 
    #     into_mach2_input_labeling $leaf_labeling $mach2_labeling

    #     # solve mach2
    #     mamba activate mach2
    #     /usr/bin/time -v \
    #         mach2 \
    #         $mach2_tree \
    #         $mach2_labeling \
    #         -o $mach2_dir \
    #         -p "LL" \
    #         -c "CM" \
    #         --max_solutions 10 \
    #         > $mach2_timing 2>&1 
    # fi

    # mamba activate metient_gpu_2
    # if [ $(ls -A "$metient_out" | wc -l) -gt 0 ]; then
    #     echo "Directory is not empty ${metient_out}"

    # echo "Processing: $data - metient"

    # if [ "$(wc -l < "$leaf_labeling")" -gt 400 ]; then
    #     echo "Too many leaves"
    #     continue
    # fi

        # /usr/bin/time -v \
        #     python \
        #     $metient \
        #     $edgelist \
        #     $leaf_labeling \
        #     -x 1 \
        #     -r "LL" \
        #     -m $metient_out \
        #     -o $metient_out \
        #     > $metient_timing 2>&1
    #     if ! compgen -G "$metient_out/*.png" > /dev/null; then
    #     #     echo $metient_out
    #         echo "No .png files found — doing something"
    #     #     # do something here




    # # convert the metient tree to a tsv
    #         awk '{print "s"$1 "\t" "s"$2}' ${metient_out}_metient_tree.txt > ${metient_out}_metient_tree.tsv
    #         echo "processing - $data - metient"

    #         python \
    #             $process_metient \
    #             $metient_out \
    #             $metient_timing \
    #             ${metient_out}_metient_tree.tsv \
    #             $edgelist \
    #             "empty"
    #     fi
    # fi

    # # metient with 3h timeout
    
    # # conda run -n metient
    # # mamba activate metient_gpu
    # # timeout -k 30s 3h \
    # #     /usr/bin/time -v \
    # #     python \
    # #     $metient \
    # #     $edgelist \
    # #     $leaf_labeling \
    # #     -o $metient_out \
    # #     -r "LL" \
    # #     > $metient_timing 2>&1 

    # # # convert the metient tree to a tsv
    # # awk '{print "s"$1 "\t" "s"$2}' ${metient_out}_metient_tree.txt > ${metient_out}_metient_tree.tsv

    # # python \
    # #     $process_metient \
    # #     $metient_out \
    # #     $metient_timing \
    # #     ${metient_out}_metient_tree.tsv \
    # #     $edgelist \
    # #     "empty"

    # solve tlp - normal
    # mamba activate tlp
    # /usr/bin/time -v \
    #     python \
    #     $tlp \
    #     fast_machina \
    #     $edgelist \
    #     $leaf_labeling \
    #     -o $tlp_normal_out \
    #     -l "LL" \
    #     > $tlp_normal_timing 2>&1 

    # solve tlp - edge regularized normal
    # mamba activate tlp
    # /usr/bin/time -v \
    #     python \
    #     $tlp \
    #     fast_machina \
    #     $edgelist \
    #     $leaf_labeling \
    #     -o $tlp_regularized_out \
    #     -e 1 \
    #     -l "LL" \
    #     > $tlp_regularized_timing 2>&1 

    # solve tlp - tree constrained
    # /usr/bin/time -v \
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