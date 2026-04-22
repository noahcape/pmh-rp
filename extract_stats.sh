#!/bin/bash
# Extract stats from the score_results.py output

score_file=$1
stats_file=$2

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
