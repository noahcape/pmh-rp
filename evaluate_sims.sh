#!/usr/bin/env bash

DIR=$1

TLP_STATS="${DIR}/tlp_stats.csv"
TLP_NORMAL_STATS="${DIR}/tlp_normal_stats.csv"
TLP_TREE_STATS="${DIR}/tlp_tree_stats.csv"
MACH2_STATS="${DIR}/mach2_stats.csv"

create_dir() {
    local new_dir=$1

    if [ ! -d $new_dir ]; then
        echo "$new_dir does not exist - making it"
        mkdir $new_dir
    fi
}

create_file() {
    local new_file=$1

    if [ ! -f $new_file ]; then
        echo "$new_file does not exist -- making it"
        touch $new_file
    fi
}

simulate_tlp_data() {
    local migration_rate=$1
    local leaves=$2
    shift 2

    local seed=$(( RANDOM % 1000 ))

    SIM_PARAMS="m${migration_rate}_l${leaves}"

    local out="${DIR}/${seed}/${SIM_PARAMS}"

    create_dir ${DIR}/${seed}

    python \
        scripts/simulations/cancer_evolution.py \
        -o $out \
        -n $leaves \
        --migration-rate $migration_rate \
        -r $seed \
        $@ # for setting an optional structure
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

extract_stats() {
    local score_file=$1
    local stats_file=$2

    edges=$(jq '.inferred_migration_graph_num_edges' $score_file)
    parsimony=$(jq '.inferred_parsimony_score' $score_file)
    tp=$(jq '.pairwise_relations.true_positives' $score_file)
    tn=$(jq '.pairwise_relations.true_negatives' $score_file)
    fp=$(jq '.pairwise_relations.false_positives' $score_file)
    fn=$(jq '.pairwise_relations.false_negatives' $score_file)
    ji=$(jq '.pairwise_relations.jaccard_index' $score_file)
    tps=$(jq '.true_parsimony_score' $score_file)
    ips=$(jq '.inferred_parsimony_score' $score_file)

    if [ ! -f $stats_file ]; then
        echo "$stats_file does not exist -- making it"
        touch $stats_file
        echo "num_edges,parsimony,tp,tn,fp,fn,ji,tps,ips" >> $stats_file
    fi

    echo "$edges,$parsimony,$tp,$tn,$fp,$fn,$ji,$tps,$ips" >> $stats_file
}

mach2_mach2sims() {
    local out=$1
    local loc=$2
    local seed=$3
    local stats_file=$4

    local out_dir="${loc}/mach2"
    local file_prefix="${loc}/T_seed${seed}"

    gtime -v \
        mach2 \
        ${file_prefix}.tree \
        ${file_prefix}.observed.labeling \
        -o $out_dir \
        --max_solutions 10 \
        > ${out}_timing.txt 2>&1 


    # loop over the results from mach2 and score each seperately
    grep '^[A-Z]-[[:space:]]*[0-9]' "${out}_timing.txt" | while read -r line; do
        first=$(echo "$line" | awk '{print $1}')
        second=$(echo "$line" | awk '{print $2}')

        inferred_labeling="${first}T-${second}.location.labeling"
        score_file=${out}_${first}${second}_scored.json

        # score result
        python \
            scripts/processing/score_result.py \
            ${file_prefix}.tree \
            ${file_prefix}.tree \
            ${file_prefix}.location.labeling \
            $out_dir/$inferred_labeling \
            ${out}_score_timing.txt \
            -o $score_file
        
        # extract stats of ancestral reconstruction
        extract_stats $score_file $stats_file

        # draw results
        dot -Tpng ${out_dir}/${first}G-${second}.dot > ${out_dir}/${first}G-${second}.png
        dot -Tpng ${out_dir}/${first}T-${second}.dot > ${out_dir}/${first}T-${second}.png
    done

    rm ${out_dir}/*.dot
    rm ${out_dir}/*.pdf

}

mach2_tlpsims() {
    local out=$1
    local loc=$2
    local seed=$3
    local stats_file=$4

    local out_dir="${loc}/mach2"
    local file_prefix="${loc}/${SIM_PARAMS}"

    gtime -v \
        mach2 \
        ${file_prefix}_perturbed.tree \
        ${file_prefix}.labeling \
        -o $out_dir \
        --max_solutions 10 \
        > ${out}_timing.txt 2>&1 


    # first count the number of solutionsfiles that mach2 will produce
    count=$(grep -E '^[0-9A-Za-z]+-[[:space:]]*[0-9]+' "${out}_timing.txt" | wc -l)


    # loop over the results from mach2 and score each seperately
    grep -E '^[0-9A-Za-z]+-[[:space:]]*[0-9]+' "${out}_timing.txt" |
    while read -r line; do
        first=$(echo "$line" | awk '{print $1}' | sed 's/-$//')   # strip trailing hyphen
        second=$(echo "$line" | awk '{print $2}')

        # prepend zero if single digit
        if [[ ${#second} -eq 1 && ${count} -gt 10 ]]; then
            second="0$second"
        fi

        inferred_labeling="${first}-T-${second}.location.labeling"
        score_file=${out}_${first}-${second}_scored.json

        # score result
        python \
            scripts/processing/score_result.py \
            ${file_prefix}_perturbed.tree \
            ${file_prefix}.tree \
            ${file_prefix}_labeling.csv \
            $out_dir/$inferred_labeling \
            ${out}_score_timing.txt \
            -o $score_file
        
        # extract stats of ancestral reconstruction
        extract_stats $score_file $stats_file

        # draw results
        dot -Tpng ${out_dir}/${first}-G-${second}.dot > ${out_dir}/${first}-G-${second}.png
        # dot -Tpng ${out_dir}/${first}-T-${second}.dot > ${out_dir}/${first}-T-${second}.png
    done

    rm ${out_dir}/*.dot
    rm ${out_dir}/*.pdf
}

tlp_mach2sims() {
    local out=$1
    local loc=$2
    local seed=$3
    local stats_file=$4

    local file_prefix="${loc}/T_seed${seed}"
    score_file=${out}_scored.json

    # solve tlp
    gtime -v \
        python \
        ./scripts/tlp.py \
        fast_machina \
        ${file_prefix}.tree \
        ${file_prefix}.observed.labeling \
        -e 1 -n 1 \
        -o $out \
        > ${out}_timing.txt 2>&1 

    # score result
    python \
        scripts/processing/score_result.py \
        ${file_prefix}.tree \
        ${file_prefix}.tree \
        ${file_prefix}.location.labeling \
        ${out}_vertex_labeling.csv \
        ${out}_score_timing.txt \
        -o $score_file

    # extract stats of ancestral reconstruction
    extract_stats $score_file $stats_file

    # draw tree
    python \
        ./scripts/plots/draw_colored_tree.py \
        ${out}_pulled_down.edgelist \
        ${out}_vertex_labeling.csv \
        -f "edgelist" \
        -o $out

    # draw results
    dot -Tpng ${out}_color_graph.dot > ${out}_color_graph.png
    dot -Tpng ${out}_colored_tree.dot > ${out}_colored_tree.png

    rm ${out}_*.dot
}

tlp_tlpsims() {
    local out=$1
    local loc=$2
    local seed=$3
    local stats_file=$4

    shift 4

    local file_prefix="${loc}/${SIM_PARAMS}"
    score_file=${out}_scored.json

    # solve tlp
    gtime -v \
        python \
        ./scripts/tlp.py \
        fast_machina \
        ${file_prefix}_perturbed_tree_edgelist.tsv \
        ${file_prefix}_leaf_labeling.csv \
        $@ \
        -o $out \
        > ${out}_timing.txt 2>&1 

    # score result
    python \
        scripts/processing/score_result.py \
        ${file_prefix}_perturbed_tree_edgelist.tsv \
        ${file_prefix}_tree_edgelist.tsv \
        ${file_prefix}_labeling.csv \
        ${out}_vertex_labeling.csv \
        ${out}_score_timing.txt \
        -o $score_file

    # extract stats of ancestral reconstruction
    extract_stats $score_file $stats_file

    # draw tree
    python \
        ./scripts/plots/draw_colored_tree.py \
        ${file_prefix}_perturbed_tree_edgelist.tsv \
        ${out}_vertex_labeling.csv \
        -f "edgelist" \
        -o $out

    # draw results
    dot -Tpng ${out}_color_graph.dot > ${out}_color_graph.png
    # dot -Tpng ${out}_colored_tree.dot > ${out}_colored_tree.png

    rm ${out}_*.dot
}

move_mach2_sims() {
    for seed in $(basename -a ${DIR}/T_seed*.location.labeling | sed -E 's/T_seed([0-9]+).*/\1/'); do
        loc=${DIR}/${seed}
        # move all the files
        mv ${DIR}/T_seed${seed}.* $loc
    done
}

produce_tlp_sims() {
    local n=$1
    local migration_rate=$2
    local leaves=$3
    shift 3

    for i in $(seq 1 "$n"); do
        simulate_tlp_data $migration_rate $leaves $@
    done
}

run_mach2_sims() {
    for loc in $DIR/*; do
        if [ -d "$loc" ]; then
            seed=$(basename "$loc")

            create_dir ${loc}

            create_dir ${loc}/tlp
            echo "Reconstructing ${seed} with tlp"
            tlp_mach2sims ${loc}/tlp/tlp_${seed} $loc $seed $TLP_STATS

            create_dir ${loc}/mach2
            echo "Reconstructing ${seed} with mach2"
            mach2_mach2sims ${loc}/mach2/mach2_${seed} $loc $seed $MACH2_STATS
        fi
    done
}

run_tlp_sims() {
    local structure=$1

    for loc in $DIR/*; do
        if [ -d $loc ]; then
            seed=$(basename $loc)

            # draw tree
            python \
                ./scripts/plots/draw_colored_tree.py \
                ${loc}/${SIM_PARAMS}_tree_edgelist.tsv \
                ${loc}/${SIM_PARAMS}_labeling.csv \
                -f "edgelist" \
                -o ${loc}/${seed}

            # draw results
            dot -Tpng ${loc}/${seed}_color_graph.dot > ${loc}/${seed}_color_graph.png
            dot -Tpng ${loc}/${seed}_colored_tree.dot > ${loc}/${seed}_colored_tree.png

            rm ${loc}/${seed}_*.dot


            true_labeling=${SIM_PARAMS}_labeling.csv
            leaf_labeling=${SIM_PARAMS}_leaf_labeling.csv
            perturbed_edgelist=${SIM_PARAMS}_perturbed_tree_edgelist.tsv
            edgelist=${SIM_PARAMS}_tree_edgelist.tsv

            echo "Reconstructing ${seed} with tlp"
            create_dir ${loc}/tlp
            tlp_tlpsims ${loc}/tlp/tlp_${seed} $loc $seed $TLP_STATS -e 1
            
            echo "Reconstructing ${seed} with tlp normal"
            create_dir ${loc}/tlp_normal
            tlp_tlpsims ${loc}/tlp_normal/tlp_normal_${seed} $loc $seed $TLP_NORMAL_STATS

            if [[ -n "$structure" ]]; then
                echo "Reconstructing ${seed} with tlp ${structure}"
                create_dir ${loc}/tlp_${structure}
                tlp_tlpsims ${loc}/tlp_${structure}/tlp_${structure}_${seed} $loc $seed $TLP_TREE_STATS -c $structure
            fi

            echo "Reconstructing ${seed} with mach2"
            create_dir ${loc}/mach2
            # turn the tlp instances into mach2 instances
            into_mach2_input_edgelist ${loc}/${perturbed_edgelist} ${loc}/${SIM_PARAMS}_perturbed.tree
            into_mach2_input_edgelist ${loc}/${edgelist} ${loc}/${SIM_PARAMS}.tree
            into_mach2_input_labeling ${loc}/${leaf_labeling} ${loc}/${SIM_PARAMS}.labeling

            mach2_tlpsims ${loc}/mach2/mach2_${seed} $loc $seed $MACH2_STATS
        fi
    done
}


# tree
produce_tlp_sims 5 0.002 1000 -s "polyclonal_tree"
run_tlp_sims "polyclonal_tree"

# dag
# produce_tlp_sims 5 0.002 1000 -s "polyclonal_dag"
# run_tlp_sims "polyclonal_dag"