# create slurm job for simulate_cancer_evolution.sh

working_dir=/n/fs/ragr-research/projects/pmh-rp
script=$working_dir/simulate_cancer_evolution.sh

job_name="simulate_cancer_evolution_tlp"
log_file=$working_dir/$job_name

sbatch -t 0-2 -c 1 -J $job_name -o $log_file.log -e $log_file.err $script