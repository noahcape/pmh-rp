#!/bin/bash
# Runs mach2 on simulated examples

seeds=(1 2 3 4 5 6 7 8 9 10)
migration_rate=(0.002)
# leaves=(50 75 100 500 750 1000)
leaves=(500 750 1000)
flips=(5 10)
# constraints=("none" "polyclonal_tree" "polyclonal_dag")
constraints=("none" "polyclonal_dag")

working_dir=/n/fs/ragr-research/projects/pmh-rp
score_labeling=${working_dir}/scripts/processing/score_result.py
out_base=${working_dir}/sims

echo "Running on $(hostname)"

extract_stats() {
    local score_file=$1
    local stats_file=$2

    if [ -f $score_file ]; then
        edges=$(jq '.inferred_migration_graph_num_edges' $score_file)
        parsimony=$(jq '.inferred_parsimony_score' $score_file)
        tp=$(jq '.pairwise_relations.true_positives' $score_file)
        tn=$(jq '.pairwise_relations.true_negatives' $score_file)
        fp=$(jq '.pairwise_relations.false_positives' $score_file)
        fn=$(jq '.pairwise_relations.false_negatives' $score_file)
        ji=$(jq '.pairwise_relations.jaccard_index' $score_file)
        tps=$(jq '.true_parsimony_score' $score_file)
        ips=$(jq '.inferred_parsimony_score' $score_file)
        time=$(jq '.elapsed_time' $score_file)

        if [ ! -f $stats_file ]; then
            echo "$stats_file does not exist -- making it"
            touch $stats_file
            echo "num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips,time" >> $stats_file
        fi

        echo "$edges,$parsimony,$tp,$tn,$fp,$fn,$ji,$tps,$ips,$time" >> $stats_file
    fi
}

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

solve_and_score() {
    local tree=$1
    # this would either be perturbed or same as above
    local solved_on_tree=$2
    local leaf_labeling=$3
    local full_labeling=$4
    local timing_file=$5
    local mach2_dir=$6
    local stats_file=$7
    local out=$8
    local score_labeling=$9

    mamba activate mach2
    # solve mach2
    /usr/bin/time -v \
        mach2 \
        $solved_on_tree \
        $leaf_labeling \
        -o $mach2_dir \
        --max_solutions 10 \
        > $timing_file 2>&1 

    # first count the number of solutionsfiles that mach2 will produce
    count=$(grep -E '^[0-9A-Za-z]+-[[:space:]]*[0-9]+' $timing_file | wc -l)

    # loop over the results from mach2 and score each seperately
    grep -E '^[0-9A-Za-z]+-[[:space:]]*[0-9]+' $timing_file |
    while read -r line; do
        first=$(echo "$line" | awk '{print $1}' | sed 's/-$//')   # strip trailing hyphen
        second=$(echo "$line" | awk '{print $2}')

        # only score the top 10 labelings
        if (( second > 10 )); then
            continue
        fi

        # prepend zero if single digit
        if [[ ${#second} -eq 1 && ${count} -gt 10 ]]; then
            second="0$second"
        fi

        inferred_labeling=${mach2_dir}/${first}-T-${second}.location.labeling
        score_file=${out}_${first}-${second}_scored.json

        mamba deactivate
        mamba activate tlp
        # score result
        python \
            $score_labeling \
            $solved_on_tree \
            $tree \
            $full_labeling \
            $inferred_labeling \
            $timing_file \
            -o $score_file
        mamba deactivate
        
        # extract stats of ancestral reconstruction
        extract_stats $score_file $stats_file
    done
}

check_dir() {
    local dir=$1
    if [ ! -d $dir ]; then
        mkdir $dir
    fi
}

export -f solve_and_score extract_stats
export score_labeling


for n_leaves in ${leaves[@]}; do
    for seed in ${seeds[@]}; do
        for constraint in ${constraints[@]}; do
            for mrate in ${migration_rate[@]}; do
                for flip in ${flips[@]}; do
                    params="leaves=${n_leaves}_constraint=${constraint}_mrate=${mrate}_flip=${flip}"
                    echo $params
                    fbase="${out_base}/${params}/$seed"
                    parent_fbase=$fbase/sim

                    # need a mach2 dir for each type of solved instance
                    mach2_dir_normal=$fbase/mach2_normal
                    check_dir $mach2_dir_normal
                    out_normal="${mach2_dir_normal}/sim_${constraint}"
                    timing_file_normal=${out_normal}_timing.txt
                    stats_file_normal="${out_base}/${params}/mach2_normal.csv"

                    mach2_dir_flip="$fbase/mach2_flip=${flip}"
                    check_dir $mach2_dir_flip
                    out_flip="${mach2_dir_flip}/sim_flip=${flip}"
                    timing_file_flip=${out_flip}_timing.txt
                    stats_file_flip="${out_base}/${params}/mach2_flip=${flip}.csv"
                    
                    mach2_dir_perturbed="$fbase/mach2_perturbed"
                    check_dir $mach2_dir_perturbed
                    out_perturbed="${mach2_dir_perturbed}/sim_perturbed"
                    timing_file_perturbed=${out_perturbed}_timing.txt
                    stats_file_perturb="${out_base}/${params}/mach2_perturbed.csv"
                    
                    leaf_labeling=${parent_fbase}_leaf_labeling.csv
                    flipped_leaf_labeling=${parent_fbase}_flipped_leaf_labeling.csv
                    perturbed_leaf_labeling=${parent_fbase}_pertubed_leaf_labeling.csv
                    full_labeling=${parent_fbase}_labeling.csv
                    perturbed_tree=${parent_fbase}_perturbed_tree_edgelist.tsv
                    tree=${parent_fbase}_tree_edgelist.tsv

                    mach2_leaf_labeling=${parent_fbase}.labeling
                    mach2_perturbed_leaf_labeling=${parent_fbase}_perturbed.labeling
                    mach2_flipped_leaf_labeling=${parent_fbase}_flipped.labeling
                    mach2_tree=${parent_fbase}.tree
                    mach2_perturbed_tree=${parent_fbase}_perturbed.tree

                    into_mach2_input_edgelist $perturbed_tree ${parent_fbase}_perturbed.tree
                    into_mach2_input_edgelist $tree ${parent_fbase}.tree
                    into_mach2_input_labeling $leaf_labeling ${parent_fbase}.labeling
                    into_mach2_input_labeling $flipped_leaf_labeling ${parent_fbase}_flipped.labeling
                    into_mach2_input_labeling $perturbed_leaf_labeling ${parent_fbase}_perturbed.labeling

                    if [[ $flip -eq 5 ]]; then
                        solve_and_score \
                            $mach2_tree \
                            $mach2_tree \
                            $mach2_leaf_labeling \
                            $full_labeling \
                            $timing_file_normal \
                            $mach2_dir_normal \
                            $stats_file_normal \
                            $out_normal \
                            $score_labeling

                        solve_and_score  \
                            $mach2_tree \
                            $mach2_perturbed_tree \
                            $mach2_perturbed_leaf_labeling \
                            $full_labeling \
                            $timing_file_perturbed \
                            $mach2_dir_perturbed \
                            $stats_file_perturb \
                            $out_perturbed \
                            $score_labeling
                    fi
                    
                    # solve_and_score \
                    #     $mach2_tree \
                    #     $mach2_tree \
                    #     $mach2_flipped_leaf_labeling \
                    #     $full_labeling \
                    #     $timing_file_flip \
                    #     $mach2_dir_flip \
                    #     $stats_file_flip \
                    #     $out_flip \
                    #     $score_labeling
                done
            done
        done
    done
done


# # solve mach2
# /usr/bin/time -v \
#     conda run -n mach2 \
#     mach2 \
#     ${parent_fbase}_perturbed.tree \
#     ${parent_fbase}.labeling \
#     -o $mach2_dir \
#     --max_solutions 10 \
#     > $timing_file 2>&1 

# # first count the number of solutionsfiles that mach2 will produce
# count=$(grep -E '^[0-9A-Za-z]+-[[:space:]]*[0-9]+' $timing_file | wc -l)

# # loop over the results from mach2 and score each seperately
# grep -E '^[0-9A-Za-z]+-[[:space:]]*[0-9]+' $timing_file |
# while read -r line; do
#     first=$(echo "$line" | awk '{print $1}' | sed 's/-$//')   # strip trailing hyphen
#     second=$(echo "$line" | awk '{print $2}')

#     # prepend zero if single digit
#     if [[ ${#second} -eq 1 && ${count} -gt 10 ]]; then
#         second="0$second"
#     fi

#     inferred_labeling=${mach2_dir}/${first}-T-${second}.location.labeling
#     score_file=${out}_${first}-${second}_scored.json

#     # score result
#     conda run -n tlp \
#         python \
#         $score_labeling \
#         ${parent_fbase}_perturbed.tree \
#         ${parent_fbase}.tree \
#         $full_labeling \
#         $inferred_labeling \
#         $timing_file \
#         -o $score_file
    
#     # extract stats of ancestral reconstruction
#     extract_stats $score_file $stats_file