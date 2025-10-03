#!/bin/bash

#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00

job_script=/home/ssoheili/genomica-code/gwas_match_with_UKB.R
module load 2024; module load R/4.4.2-gfbf-2024a

Rscript $job_script

