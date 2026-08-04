#!/bin/bash

#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=10:00:00

module load 2024; module load R/4.4.2-gfbf-2024a
Rscript /home/ssoheili/genomica-code/psy_match_with_UKB.R 
