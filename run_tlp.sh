#!/bin/bash
# Runs tree-labeling-polytope on simulated examples
set -euo pipefail
IFS=$'\n\t'

seeds=(1 2 3 4 5 6 7 8 9 10)
migration_rate=(0.002)
# leaves=(50 75 100 500 750 1000)
leaves=(500 750 1000)
flips=(5 10)
constraints=("polyclonal_tree")
# constraints=("none" "polyclonal_dag")
# constraints=("none" "polyclonal_tree" "polyclonal_dag")

working_dir=/n/fs/ragr-research/projects/pmh-rp
tlp=${working_dir}/scripts/tlp.py
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

solve_and_score() {
    local tree=$1
    local solved_on_tree=$2
    local leaf_labeling=$3
    local full_labeling=$4
    local regularize=$5
    local constraint=$6
    local timing_file=$7
    local score_file=$8
    local stats_file=$9
    local out="${10}"

    # solve tlp
    /usr/bin/time -v \
        python \
        $tlp \
        fast_machina \
        $solved_on_tree \
        $leaf_labeling \
        -o $out \
        -e $regularize \
        -c $constraint \
        > $timing_file 2>&1 

    inferred_labeling=${out}_vertex_labeling.csv

    # score result
    python \
        $score_labeling \
        $solved_on_tree \
        $tree \
        $full_labeling \
        $inferred_labeling \
        $timing_file \
        -o $score_file

    # extract stats of ancestral reconstruction
    extract_stats $score_file $stats_file
}

mamba activate tlp

# TODO: run normal tlp labeling on everything
for n_leaves in ${leaves[@]}; do
    for seed in ${seeds[@]}; do
        for constraint in ${constraints[@]}; do
            for mrate in ${migration_rate[@]}; do
                for flip in ${flips[@]}; do
                    params="leaves=${n_leaves}_constraint=${constraint}_mrate=${mrate}_flip=${flip}"
                    echo "solving $params"

                    fbase="${out_base}/${params}/$seed/sim"
                    leaf_labeling=${fbase}_leaf_labeling.csv
                    flipped_leaf_labeling=${fbase}_flipped_leaf_labeling.csv
                    perturbed_leaf_labeling=${fbase}_pertubed_leaf_labeling.csv
                    full_labeling=${fbase}_labeling.csv
                    perturbed_tree=${fbase}_perturbed_tree_edgelist.tsv
                    tree=${fbase}_tree_edgelist.tsv


                    # for each run on normal, flipped, and perturbed
                    # if constraint is none -> normal, regularized
                    # if constraint is dag, tree -> normal, regularized, constrained
                    # only do normal and perturbed when flip = 5

                    # flipped
                    out_normal_flip="${fbase}_tlp_normal_flip=${flip}"
                    timing_file_normal_flip=${out_normal_flip}_timing.txt                    
                    score_file_normal_flip=${out_normal_flip}_results.json
                    stats_file_normal_flip="${out_base}/${params}/tlp_normal_flip=${flip}.csv"

                    solve_and_score \
                        $tree \
                        $tree \
                        $flipped_leaf_labeling \
                        $full_labeling \
                        0 \
                        "none" \
                        $timing_file_normal_flip \
                        $score_file_normal_flip \
                        $stats_file_normal_flip \
                        $out_normal_flip

                    out_reg_flip="${fbase}_tlp_reg_flip=${flip}"
                    timing_file_reg_flip="${out_reg_flip}_timing.txt"
                    score_file_reg_flip="${out_reg_flip}_results.json"
                    stats_file_reg_flip="${out_base}/${params}/tlp_reg_flip=${flip}.csv"

                    solve_and_score \
                        $tree \
                        $tree \
                        $flipped_leaf_labeling \
                        $full_labeling \
                        1 \
                        "none" \
                        $timing_file_reg_flip \
                        $score_file_reg_flip \
                        $stats_file_reg_flip \
                        $out_reg_flip

                    if [[ "$constraint" != "none" ]]; then
                        out_constrained_flip="${fbase}_tlp_${constraint}_flip=${flip}"
                        timing_file_constrained_flip="${out_constrained_flip}_timing.txt"
                        score_file_constrained_flip="${out_constrained_flip}_results.json"
                        stats_file_constrained_flip="${out_base}/${params}/tlp_${constraint}_flip=${flip}.csv"

                        solve_and_score \
                            $tree \
                            $tree \
                            $flipped_leaf_labeling \
                            $full_labeling \
                            0 \
                            $constraint \
                            $timing_file_constrained_flip \
                            $score_file_constrained_flip \
                            $stats_file_constrained_flip \
                            $out_constrained_flip
                    fi

                    if [[ $flip -eq 5 ]]; then
                        # normal
                        out_normal="${fbase}_tlp_normal"
                        timing_file_normal=${out_normal}_timing.txt                    
                        score_file_normal=${out_normal}_results.json
                        stats_file_normal="${out_base}/${params}/tlp_normal.csv"

                        solve_and_score \
                            $tree \
                            $tree \
                            $leaf_labeling \
                            $full_labeling \
                            0 \
                            "none" \
                            $timing_file_normal \
                            $score_file_normal \
                            $stats_file_normal \
                            $out_normal

                        out_reg="${fbase}_tlp_reg"
                        timing_file_reg="${out_reg}_timing.txt"
                        score_file_reg="${out_reg}_results.json"
                        stats_file_reg="${out_base}/${params}/tlp_reg.csv"

                        solve_and_score \
                            $tree \
                            $tree \
                            $leaf_labeling \
                            $full_labeling \
                            1 \
                            "none" \
                            $timing_file_reg \
                            $score_file_reg \
                            $stats_file_reg \
                            $out_reg

                        if [[ "$constraint" != "none" ]]; then
                            out_constrained="${fbase}_tlp_${constraint}"
                            timing_file_constrained="${out_constrained}_timing.txt"
                            score_file_constrained="${out_constrained}_results.json"
                            stats_file_constrained="${out_base}/${params}/tlp_${constraint}.csv"

                            solve_and_score \
                                $tree \
                                $tree \
                                $leaf_labeling \
                                $full_labeling \
                                0 \
                                $constraint \
                                $timing_file_constrained \
                                $score_file_constrained \
                                $stats_file_constrained \
                                $out_constrained
                        fi

                        # perturbed
                        out_normal_perturbed="${fbase}_tlp_normal_perturbed"
                        timing_file_normal_perturbed=${out_normal_perturbed}_timing.txt                    
                        score_file_normal_perturbed=${out_normal_perturbed}_results.json
                        stats_file_normal_perturbed="${out_base}/${params}/tlp_normal_perturbed.csv"

                        solve_and_score \
                            $tree \
                            $perturbed_tree \
                            $perturbed_leaf_labeling \
                            $full_labeling \
                            0 \
                            "none" \
                            $timing_file_normal_perturbed \
                            $score_file_normal_perturbed \
                            $stats_file_normal_perturbed \
                            $out_normal_perturbed

                        out_reg_perturbed="${fbase}_tlp_reg_perturbed"
                        timing_file_reg_perturbed="${out_reg_perturbed}_timing.txt"
                        score_file_reg_perturbed="${out_reg_perturbed}_results.json"
                        stats_file_reg_perturbed="${out_base}/${params}/tlp_reg_perturbed.csv"

                        solve_and_score \
                            $tree \
                            $perturbed_tree \
                            $perturbed_leaf_labeling \
                            $full_labeling \
                            1 \
                            "none" \
                            $timing_file_reg_perturbed \
                            $score_file_reg_perturbed \
                            $stats_file_reg_perturbed \
                            $out_reg_perturbed

                        if [[ "$constraint" != "none" ]]; then
                            out_constrained_perturbed="${fbase}_tlp_${constraint}_perturbed"
                            timing_file_constrained_perturbed="${out_constrained_perturbed}_timing.txt"
                            score_file_constrained_perturbed="${out_constrained_perturbed}_results.json"
                            stats_file_constrained_perturbed="${out_base}/${params}/tlp_${constraint}_perturbed.csv"

                            solve_and_score \
                                $tree \
                                $perturbed_tree \
                                $perturbed_leaf_labeling \
                                $full_labeling \
                                0 \
                                $constraint \
                                $timing_file_constrained_perturbed \
                                $score_file_constrained_perturbed \
                                $stats_file_constrained_perturbed \
                                $out_constrained_perturbed
                        fi
                    fi                    
                done
            done
        done
    done
done


                    # solve tlp
                    # /usr/bin/time -v \
                    #     conda run -n tlp \
                    #     python \
                    #     $tlp \
                    #     fast_machina \
                    #     $pertubed_tree \
                    #     $leaf_labeling \
                    #     -o $out \
                    #     -e $regularize \
                    #     > $timing_file 2>&1 
                    #     # -c $constraint \

                    # inferred_labeling=${out}_vertex_labeling.csv

                    # # score result
                    # conda run -n tlp \
                    #     python \
                    #     $score_labeling \
                    #     $pertubed_tree \
                    #     $tree \
                    #     $full_labeling \
                    #     $inferred_labeling \
                    #     $timing_file \
                    #     -o $score_file

                    # # extract stats of ancestral reconstruction
                    # extract_stats $score_file $stats_file