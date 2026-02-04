#!/bin/bash

#SBATCH --partition=fat_genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1440G
#SBATCH --time=72:00:00

module load 2024
module load python-isal/1.7.0-GCCcore-13.3.0
#python -u h5_concat.py
python -u h5_mask_concat.py
