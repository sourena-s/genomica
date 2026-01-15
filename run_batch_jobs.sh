#!/bin/bash

#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=120:00:00


n_total=108        # total number of array jobs
n_batch=2         # number of array jobs per batch (each 16 cores)
script="panukb_dl_slurm.sh"

for (( array_id=0; array_id<n_total; array_id+=1 )); do
    while true; do
	    if [ "$(squeue -u ssoheili |grep panukb| wc -l)" -lt $n_batch ]; then
            break
        fi
        echo "Waiting for jobs to finish... $(squeue -u ssoheili |grep panukb| wc -l) jobs in queue)"
        sleep 120
    done

    echo "Submitting batch $array_id"
    sbatch -a "$array_id" "$script"
    sleep 120

done

