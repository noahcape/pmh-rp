#!/usr/bin/env bash

DIR=$1
OUT=$2
LEAVES=${3:-1000}

simulate_data() {
    local out=$1
    local mig_rate=$2
    local random_num=$3
    local leaves=$4
    local structure=$5

    echo "Simulation data"
    python \
        scripts/simulations/cancer_evolution.py \
        -o $out \
        -n $leaves \
        --migration-rate $mig_rate \
        -r $random_num \
        -s $structure
}

solve_tlp() {
    local mode=$1
    local outfile=$2
    local dir=$3

    shift 3

    echo "Solving TLP - $mode"
    gtime -v \
        python \
        scripts/tlp.py \
        fast_machina \
        $dir/${OUT}_perturbed_tree_edgelist.tsv \
        $dir/${OUT}_leaf_labeling.csv \
        $@ \
        -o $dir/$outfile \
        > ${dir}/timing_${outfile}.txt 2>&1
}

solve_tlp_no_perturb() {
    local mode=$1
    local outfile=$2
    local dir=$3

    shift 3

    echo "Solving TLP - $mode"
    python \
        scripts/tlp.py \
        fast_machina \
        $dir/${OUT}_tree_edgelist.tsv \
        $dir/${OUT}_leaf_labeling.csv \
        $@ \
        -o $dir/$outfile
}

score_result() {
    local mode=$1
    local outfile=$2
    local score_json=$3
    local dir=$4

    echo "Scoring results - $mode"
    python \
        scripts/processing/score_result.py \
        $dir/${OUT}_perturbed_tree_edgelist.tsv \
        $dir/${OUT}_tree_edgelist.tsv \
        $dir/${OUT}_labeling.csv \
        $dir/${outfile}_vertex_labeling.csv \
        $dir/timing.txt \
        -o $score_json
}

score_result_no_perturb() {
    local mode=$1
    local outfile=$2
    local score_json=$3
    local dir=$4

    echo "Scoring results - $mode"
    python \
        scripts/processing/score_result.py \
        $dir/${OUT}_tree_edgelist.tsv \
        $dir/${OUT}_tree_edgelist.tsv \
        $dir/${OUT}_labeling.csv \
        $dir/${outfile}_vertex_labeling.csv \
        $dir/timing.txt \
        -o $score_json
}

draw_tree() {
    local mode=$1
    local outfile=$2
    local labeling=$3
    local dir=$4

    echo "Drawing tree -- $mode"
    python \
        scripts/plots/draw_colored_tree.py \
        $dir/${OUT}_tree_edgelist.tsv \
        $labeling \
        -f "edgelist" \
        -o $outfile

    dot -Tpng ${outfile}_color_graph.dot \
        > ${outfile}_color_graph.png
    dot -Tpng ${outfile}_colored_tree.dot \
        > ${outfile}_colored_tree.png
}

extract_stats() {
    local outfile=$1
    local allstats=$2

    echo "Extracting relevant stats"
    local edges
    local parsimony
    local fpr

    edges=$(jq '.inferred_migration_graph_num_edges' $outfile)
    parsimony=$(jq '.inferred_parsimony_score' $outfile)
    fpr=$(jq '.pairwise_relations.false_positive_rate' $outfile)

    if [ ! -f $allstats ]; then
        echo "$allstats does not exist -- making it"
        touch $allstats
        echo "num_edges,parsimony,fpr" >> $allstats
    fi

    echo "$edges,$parsimony,$fpr" >> $allstats
}

# solve regularized parsimony on one large instance with different
# regularized weighting value \lambda
find_pareto_front() {
    # reg_weights=(0 0.0001 0.001 0.01 0.1 0.5 1 2 3)
    reg_weights=(0.3 0.8)

    seed=42
    dir="${DIR}/automate_regularization_param"

    if [ ! -d $dir ]; then
        echo "$dir does not exist - making it"
        mkdir $dir
    fi

    #simulate one data set
    simulate_data $dir/$OUT 0.0015 $seed $LEAVES none
    draw_tree "simulated" $dir/${out_sim} $dir/${OUT}_labeling.csv $dir


    out_reg=inferred_labeling_regularized
    out_sim="ground_truth_sim"
    out_with_reg_json=${dir}/${OUT}_results_with_reg.json

    # $dir/$out_reg
    # solve tree labeling polytope with regularization and some regularization weight
    # this will write a json file which contains the object, migrations, edges, and regularization weight
    solve_tlp "regularized parsimony" $out_reg $dir -e 1

    # score the inferred labeling against the simulated tree
    score_result "with regularization" $out_reg $out_with_reg_json $dir

    # draw both the inferred tree through regularization as well as the simulated tree to compare
    draw_tree "with regularization" $dir/${out_reg} $dir/${out_reg}_vertex_labeling.csv $dir
}

simulate_many() {
    local structure=$1

    for i in {1..5}; do
        num=$(( 10 + RANDOM % 100 ))

        local dir=${DIR}_${i}_${structure}

        if [ ! -d $dir ]; then
            echo "$dir does not exist - making it"
            mkdir $dir
        fi

        # modify the rate
        simulate_data $dir/$OUT 0.008 $num $LEAVES $structure

        touch $dir/timing.txt

        out_normal="inferred_labeling_normal"
        out_regularized="inferred_labeling_regularized"
        out_additional_label="inferred_labeling_additional_label_regularized"
        out_sim="ground_truth_sim"

        # structure choices "polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag", "none"
        out_polyclonal_tree="inferred_labeling_polyclonal_tree"
        out_polyclonal_dag="inferred_labeling_polyclonal_dag"

        out_no_reg_json=${dir}/${OUT}_results_no_reg.json
        out_with_reg_json=${dir}/${OUT}_results_with_reg.json
        
        out_polyclonal_tree_json="${dir}/${OUT}_results_polyclonal_tree.json"
        out_polyclonal_dag_json="${dir}/${OUT}_results_polyclonal_dag.json"


        # for the given structure solve tlp as fastMACHINA
        solve_tlp "normal parsimony" $out_normal $dir
        solve_tlp "regularized parsimony" $out_regularized $dir -e 1

        if [[ "$structure" == "polyclonal_tree" ]]; then
            solve_tlp $structure $out_polyclonal_tree $dir -c ${structure}
            score_result $structure $out_polyclonal_tree $out_polyclonal_tree_json $dir
            draw_tree $structure "$dir/$out_polyclonal_tree" "$dir/${out_polyclonal_tree}_vertex_labeling.csv" $dir
        fi

        if [[ "$structure" == "polyclonal_dag" ]]; then
            solve_tlp $structure $out_polyclonal_dag $dir -c ${structure}
            score_result $structure $out_polyclonal_dag $out_polyclonal_dag_json $dir
            draw_tree $structure "$dir/$out_polyclonal_dag" "$dir/${out_polyclonal_dag}_vertex_labeling.csv" $dir
        fi
        
        score_result "without regularization" $out_normal $out_no_reg_json $dir
        score_result "with regularization" $out_regularized $out_with_reg_json $dir

        draw_tree "no regularization" $dir/${out_normal} $dir/${out_normal}_vertex_labeling.csv $dir
        draw_tree "with regularization" $dir/${out_regularized} $dir/${out_regularized}_vertex_labeling.csv $dir
        draw_tree "simulated" $dir/${out_sim} $dir/${OUT}_labeling.csv $dir

        # extract_stats $out_no_reg_json "$DIR/results_no_reg_weight_08.csv"
        # extract_stats $out_with_reg_json "$DIR/results_with_reg.csv"

        echo "Clear dot files"
        rm $dir/*.dot
    done
}

execute_one() {
    local labeling=$1
    local edge_list=$2

    local out_normal="inferred_labeling_normal"
    local out_regularized="inferred_labeling_regularized"
    local out_pdag="inferred_labeling_polyclonal_dag"
    local out_ptree="inferred_labeling_polyclonal_tree"

    solve_tlp_no_perturb "normal parsimony" $out_normal $DIR -l LL
    solve_tlp_no_perturb "regularized parsimony" $out_regularized $DIR -e 1 -l LL
    solve_tlp_no_perturb "polyclonal dag" $out_pdag $DIR -l LL -c polyclonal_dag
    solve_tlp_no_perturb "polyclonal tree" $out_ptree $DIR -l LL -c polyclonal_tree

    draw_tree "no regularization" $DIR/${out_normal} $DIR/${out_normal}_vertex_labeling.csv $DIR
    draw_tree "with regularization" $DIR/${out_regularized} $DIR/${out_regularized}_vertex_labeling.csv $DIR
    draw_tree "polyclonal dag" $DIR/${out_pdag} $DIR/${out_pdag}_vertex_labeling.csv $DIR
    draw_tree "with ptree" $DIR/${out_ptree} $DIR/${out_ptree}_vertex_labeling.csv $DIR
}

collect_statistics() {
    local structure=$1
    for i in {1..10}; do
        file_normal=examples/cancer_evolution/simulations/results_${i}_${structure}/sim_results_${structure}.json

        edges=$(jq '.inferred_migration_graph_num_edges' $file_normal)
        parsimony=$(jq '.inferred_parsimony_score' $file_normal)
        tp=$(jq '.pairwise_relations.true_positives' $file_normal)
        tn=$(jq '.pairwise_relations.true_negatives' $file_normal)
        fp=$(jq '.pairwise_relations.false_positives' $file_normal)
        fn=$(jq '.pairwise_relations.false_negatives' $file_normal)
        ji=$(jq '.pairwise_relations.jaccard_index' $file_normal)

        tps=$(jq '.true_parsimony_score' $file_normal)
        ips=$(jq '.inferred_parsimony_score' $file_normal)

        allstats_normal=examples/cancer_evolution/simulations/simulations_${structure}_results_${structure}.csv
        if [ ! -f $allstats_normal ]; then
            echo "$allstats_normal does not exist -- making it"
            touch $allstats_normal
            echo "num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips" >> $allstats_normal
        fi

        echo "$edges,$parsimony,$tp,$tn,$fp,$fn,$ji,$tps,$ips" >> $allstats_normal


        file_regularized=examples/cancer_evolution/simulations/results_${i}_${structure}/sim_results_with_reg.json

        edges=$(jq '.inferred_migration_graph_num_edges' $file_regularized)
        parsimony=$(jq '.inferred_parsimony_score' $file_regularized)
        fp=$(jq '.pairwise_relations.false_positives' $file_regularized)
        fn=$(jq '.pairwise_relations.false_negatives' $file_regularized)
        tn=$(jq '.pairwise_relations.true_negatives' $file_regularized)
        tp=$(jq '.pairwise_relations.true_positives' $file_regularized)
        ji=$(jq '.pairwise_relations.jaccard_index' $file_regularized)

        tps=$(jq '.true_parsimony_score' $file_regularized)
        ips=$(jq '.inferred_parsimony_score' $file_regularized)

        allstats_regularized=examples/cancer_evolution/simulations/simulations_${structure}_results_large_regularized.csv
        if [ ! -f $allstats_regularized ]; then
            echo "$allstats_regularized does not exist -- making it"
            touch $allstats_regularized
            echo "num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips" >> $allstats_regularized
        fi

        echo "$edges,$parsimony,$tp,$tn,$fp,$fn,$ji,$tps,$ips" >> $allstats_regularized

        file_base=examples/cancer_evolution/simulations/results_${i}_${structure}/sim_results_no_reg.json

        edges=$(jq '.inferred_migration_graph_num_edges' $file_base)
        parsimony=$(jq '.inferred_parsimony_score' $file_base)
        fp=$(jq '.pairwise_relations.false_positives' $file_base)
        fn=$(jq '.pairwise_relations.false_negatives' $file_base)
        tn=$(jq '.pairwise_relations.true_negatives' $file_base)
        tp=$(jq '.pairwise_relations.true_positives' $file_base)
        ji=$(jq '.pairwise_relations.jaccard_index' $file_base)

        tps=$(jq '.true_parsimony_score' $file_base)
        ips=$(jq '.inferred_parsimony_score' $file_base)

        allstats_file_base=examples/cancer_evolution/simulations/simulations_${structure}_results_large_no_reg.csv
        if [ ! -f $allstats_file_base ]; then
            echo "$allstats_file_base does not exist -- making it"
            touch $allstats_file_base
            echo "num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips" >> $allstats_file_base
        fi

        echo "$edges,$parsimony,$tp,$tn,$fp,$fn,$ji,$tps,$ips" >> $allstats_file_base

    done
}

# execute a single instance -- here we reanalyze the CP28 clone comparing normal to regularized
# execute_one "./examples/cancer_evolution/CP28_leaf_labeling.csv" "./examples/cancer_evolution/CP28_tree_edgelist.tsv"

# simulate cancer evolution then compare normal parsimony to regularized parsimony
# # structure choices "polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag", "none"
# simulate_many
# simulate_many "polyclonal_tree"
# simulate_many "polyclonal_dag"

# solve regularized parsimony 10 times with different regularization weights
find_pareto_front

# collect_statistics "polyclonal_tree"
# collect_statistics "polyclonal_dag"

