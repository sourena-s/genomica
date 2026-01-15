#!/bin/bash

#SBATCH --partition=fat_genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1400G
#SBATCH --time=72:00:00


export FSLOUTPUTTYPE=NIFTI2_GZ
melodic=/projects/0/genica/sourena/software/fsl/bin/melodic

dim=4

input_file=/home/ssoheili/genomica-code/matrix_4d_psyr2.nii.gz

input_file=/home/ssoheili/genomica-code/matrix_4d_manual.nii.gz
output_path=${input_file%.nii.gz}_ICA_dim_$dim
 
mkdir -p $output_path

$melodic -i $input_file --update_mask -o $output_path --Oall --Opca -d $dim --nobet --report -v --keep_meanvol --vn

