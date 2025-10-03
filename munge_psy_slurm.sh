#!/bin/bash

#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
######SBATCH --exclusive
#SBATCH --time=10:00:00

job_script=/home/ssoheili/genomica-code/munge_psy.sh

bash $job_script
