#!/bin/bash -l

#SBATCH -A wren-lab
#SBATCH --output slurm-%x.%A.%a.log
##SBATCH --mail-user pandas@omrf.org
#SBATCH --mail-type END,FAIL,ARRAY_TASKS
#SBATCH --job-name genom-runner-lite
#SBATCH --mem 4G
#SBATCH -p serial
#SBATCH --nodes 1
#SBATCH --cpus-per-task 10
#SBATCH -t 0:1:15
module load python/3.12.3
python3 main.py
