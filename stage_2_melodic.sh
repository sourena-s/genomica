#!/bin/bash
#SBATCH --partition=fat_genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1440G
#SBATCH --time=72:00:00


##FSLDIR=/projects/0/genica/sourena/software/fsl
##PATH=${FSLDIR}/share/fsl/bin:${PATH}
##export FSLDIR PATH
##. ${FSLDIR}/etc/fslconf/fsl.sh
module load 2025
module load Anaconda3/2025.06-1
module load OpenBLAS/0.3.29-GCC-14.2.0
conda activate fsl-melodic
export FSLOUTPUTTYPE=NIFTI2_GZ
#export FSLOUTPUTTYPE=NIFTI_GZ
#melodic=/projects/0/genica/sourena/software/fsl/bin/melodic
melodic=/home/ssoheili/.conda/envs/fsl-melodic/bin/melodic

dim=50

input_file=/home/ssoheili/genomica-code/matrix_4d_r2psy.nii.gz
#input_file=/home/ssoheili/genomica-code/matrix_4d_manual.nii.gz
#input_file=/home/ssoheili/genomica-code/matrix_4d_manual_rpoint5.nii.gz
#input_file=/home/ssoheili/genomica-code/matrix_4d_manual_rpoint3.nii.gz
#input_file=/home/ssoheili/genomica-code/matrix_4d_manual_rpoint3_maf_1percent.nii.gz

output_path=${input_file%.nii.gz}_ICA_dim_$dim
 
mkdir -p $output_path

echo "Free memory in GB: $(free -g)"

$melodic -i $input_file --update_mask -o $output_path --Oall --Opca -d $dim --nobet --no_mm -v --keep_meanvol --vn

