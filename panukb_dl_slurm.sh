#!/bin/bash

#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
######SBATCH --exclusive
#SBATCH --time=10:00:00

job_script=/home/ssoheili/genomica-code/panukb_dl.sh

bash $job_script $SLURM_ARRAY_TASK_ID 16


####after running all jobs concatenate their rho and z
## cat panukbb/sumstats/EUR/*psy_rg.rho|awk 'NF==13 && (NR==1 || ($1!=""  && $1!="Phenotype"))' > panukbb/sumstats/pheno_psy_gencor.rho
## cat panukbb/sumstats/EUR/*psy_rg.z  |awk 'NF==13 && (NR==1 || ($1!=""  && $1!="Phenotype"))' > panukbb/sumstats/pheno_psy_gencor.z




