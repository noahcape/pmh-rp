#!/bin/bash

# Simulates cancer evolution using ./scripts/simulations/cancer_evolution.py

seeds=(1 2 3 4 5 6 7 8 9 10)
migration_rate=(0.001)
leaves=(50 100 500 750 1000)
# flips=(5 10)
# constraints=("none" "polyclonal_tree" "polyclonal_dag")
constraints=("none")

working_dir=/n/fs/ragr-research/projects/pmh-rp
simulation_file=${working_dir}/scripts/simulations/cancer_evolution.py
out_base=${working_dir}/new_sims

for seed in ${seeds[@]}; do
    for n_leaves in ${leaves[@]}; do
        for constraint in ${constraints[@]}; do
            for mrate in ${migration_rate[@]}; do
                params="leaves=${n_leaves}_mrate=${mrate}"
                
                if [ ! -d ${out_base}/${params} ]; then
                    mkdir ${out_base}/${params}
                fi

                if [ ! -d ${out_base}/${params}/${seed} ]; then
                    mkdir ${out_base}/${params}/${seed}
                fi

                fout="${out_base}/${params}/$seed/sim"

                conda run -n tlp \
                    python $simulation_file \
                    -o $fout \
                    -n $n_leaves \
                    -r $seed \
                    -s "none" \
                    --migration-rate $mrate
            done
        done
    done
done
