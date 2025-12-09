#!/bin/bash
#SBATCH --job-name=star_alignment_pipe
#SBATCH -p cpu
#SBATCH --nodes=1
# SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -o logs/star_alignment_pipe.%N.%j.out 
#SBATCH -e logs/star_alignment_pipe.%N.%j.err 
#SBATCH --array=1-8

arrayjob=`cat star_alignment_pipe.sh | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

bash -c "$arrayjob && {     echo 'Finish successfully! '; } || {     echo $arrayjob >> star_alignment_pipe.failed.sh; }"
echo End time : `date`
