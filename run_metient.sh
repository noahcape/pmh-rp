#!/bin/bash
# Runs metient on simulated examples

seeds=(1 2 3 4 5 6 7 8 9 10)
# seeds=(7 8 9 10)
migration_rate=(0.002)
# only run on 50 75 100
leaves=(500)
#  100)
#  500 750 1000)
flips=(5)
# constraints=("none" "polyclonal_tree" "polyclonal_dag")
constraints=("none" "polyclonal_dag")

working_dir=/n/fs/ragr-research/projects/pmh-rp
metient=${working_dir}/scripts/processing/run_metient.py
process_metient=${working_dir}/scripts/processing/process_metient_output.py
out_base=${working_dir}/sims

echo "Running on $(hostname)"

solve_and_process() {
    local metient=$1
    local tree=$2
    local solved_on_tree=$3
    local leaf_labeling=$4
    local full_labeling=$5
    local metient_out=$6
    local seed=$7
    local timing_file=$8
    local process_metient=$9
    local stats_file="${10}"
    local fbase="${11}"

    # solve metient use a 5 hour timelimit
        # conda run -n metient \
    /usr/bin/time -v \
        python \
        $metient \
        $solved_on_tree \
        $leaf_labeling \
        -f $full_labeling \
        -o $fbase \
        -x $seed \
        -m $metient_out \
        > $timing_file 2>&1 

    # put the score file in the specific dir for the trial (flipped, normal, perturbed)
    # put the stats file in the base directory of the parameters
    # process output

    # parser.add_argument("metient", help="base directory to find files") # this is going to be the metient dir
    # parser.add_argument("timing", help="Timing file from running metient")
    # parser.add_argument("tree", help="Ground truth tree")
    # parser.add_argument("solved_tree", help="Tree solved reconstruction on")
    # parser.add_argument("labeling", help="Full labeling of tree")
    # parser.add_argument("stats", help="File to collect stats")
    # conda run -n metient \
    python \
        $process_metient \
        $metient_out \
        $timing_file \
        $tree \
        $solved_on_tree \
        $full_labeling \
        -s $stats_file
}

export -f solve_and_process
export metient
export process_metient


mamba activate metient_gpu
for n_leaves in ${leaves[@]}; do
    for seed in ${seeds[@]}; do
        for constraint in ${constraints[@]}; do
            for mrate in ${migration_rate[@]}; do
                for flip in ${flips[@]}; do
                    params="leaves=${n_leaves}_constraint=${constraint}_mrate=${mrate}_flip=${flip}"
                    echo "$params"

                    # data for solving and processing output 
                    fbase="${out_base}/${params}/$seed/sim"
                    leaf_labeling=${fbase}_leaf_labeling.csv
                    flipped_leaf_labeling=${fbase}_flipped_leaf_labeling.csv
                    perturbed_leaf_labeling=${fbase}_pertubed_leaf_labeling.csv
                    full_labeling=${fbase}_labeling.csv
                    perturbed_tree=${fbase}_perturbed_tree_edgelist.tsv
                    tree=${fbase}_tree_edgelist.tsv

                    # flipped
                    # metient_flipped="${out_base}/${params}/${seed}/metient_flip=${flip}"
                    # timing_file_flipped="${fbase}_metient_flip=${flip}_timing.txt"
                    # stats_file_flipped="${out_base}/${params}/metient_flip=${flip}.csv"

                    # timeout --foreground --signal=SIGTERM --kill-after=5m 24h \
                    #     bash -c 'solve_and_process "$@"' _ \
                    #     $metient \
                    #     $tree \
                    #     $tree \
                    #     $flipped_leaf_labeling \
                    #     $full_labeling \
                    #     $metient_flipped \
                    #     $seed \
                    #     $timing_file_flipped \
                    #     $process_metient \
                    #     $stats_file_flipped \
                    #     $fbase \
                    #     || echo "Timed out (expected)"

                    if [[ $flip -eq 5 ]]; then
                        # normal
                        metient_normal="${out_base}/${params}/${seed}/metient_normal"
                        timing_file_normal="${fbase}_metient_normal_timing.txt"
                        stats_file_normal="${out_base}/${params}/metient_normal.csv"

                        timeout --foreground --signal=SIGTERM --kill-after=5m 24h \
                            bash -c 'solve_and_process "$@"' _ \
                            $metient \
                            $tree \
                            $tree \
                            $leaf_labeling \
                            $full_labeling \
                            $metient_normal \
                            $seed \
                            $timing_file_normal \
                            $process_metient \
                            $stats_file_normal \
                            $fbase \
                            || echo "Timed out (expected)"

                        # perturbed
                        metient_perturbed="${out_base}/${params}/${seed}/metient_perturbed"
                        timing_file_perturbed="${fbase}_metient_perturbed_timing.txt"
                        stats_file_perturbed="${out_base}/${params}/metient_perturbed.csv"

                        timeout --foreground --signal=SIGTERM --kill-after=5m 24h \
                            bash -c 'solve_and_process "$@"' _ \
                            $metient \
                            $tree \
                            $perturbed_tree \
                            $perturbed_leaf_labeling \
                            $full_labeling \
                            $metient_perturbed \
                            $seed \
                            $timing_file_perturbed \
                            $process_metient \
                            $stats_file_perturbed \
                            $fbase \
                            || echo "Timed out (expected)"
                    fi
                done
            done
        done
    done
done