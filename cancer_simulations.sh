sim="/n/fs/ragr-research/projects/cancer_sims/target/release/cancer_migration_sims"
sim_dir=rust_sims_v2
seeds=(1 2 3 4 5 6 7 8 9 10)
#  6 7 8 9 10)
g=(6 7 8 9 10 11)
migration_rate=0.01

if [[ ! -d $sim_dir ]]; then 
    mkdir $sim_dir
fi

for gens in ${g[@]}; do
    for seed in ${seeds[@]}; do
        out="${sim_dir}/${gens}/${seed}/sim"
        if [[ ! -d "${sim_dir}/${gens}" ]]; then 
            mkdir "${sim_dir}/${gens}"
        fi

        if [[ ! -d "${sim_dir}/${gens}/${seed}" ]]; then 
            mkdir "${sim_dir}/${gens}/${seed}"
        fi
        
        $sim \
            -o $out \
            -g $gens \
            -s 6 \
            -m $migration_rate \
            -r $seed \
            -t 0.25 \
            -p 3 \
            > /dev/null
    done
done