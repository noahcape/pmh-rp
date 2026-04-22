# create slurm job for run_tlp.sh

working_dir=/n/fs/ragr-research/projects/pmh-rp
script=$working_dir/run_tlp.sh

job_name="solve_tlp"
log_file=$working_dir/$job_name

sbatch -t 1-0 -c 1 -J $job_name -o $log_file.log -e $log_file.err $script


