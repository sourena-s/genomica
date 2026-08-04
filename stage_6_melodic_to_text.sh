#!/bin/bash
module load 2025; module load Anaconda3/2025.06-1
conda activate fslavw

catalogue=/home/ssoheili/genomica-code/clumps/merged_manual-r2point3_psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue
genomica1=/projects/0/einf2700/sourena/genica-2025/melodic_results/matrix_4d_manual_mergedmask_mergedmask_ICA_dim_20_nl_pow3/melodic_oIC.nii.gz
genomica2=/projects/0/einf2700/sourena/genica-2025/melodic_results/matrix_4d_r2psy_mergedmask_mergedmask_ICA_dim_20_nl_pow3/melodic_oIC.nii.gz


dim=20

fsl2ascii $genomica1 ${genomica1%.nii.gz}.txt
fsl2ascii $genomica2 ${genomica2%.nii.gz}.txt

for ((i=0;i<20;i++));do
	cat "${genomica1%.nii.gz}.txt$(printf %05d $i)"|tr " " "\n" |awk '$1!=""'  > "${genomica1%.nii.gz}.ic$((i+1))"
	cat "${genomica2%.nii.gz}.txt$(printf %05d $i)"|tr " " "\n" |awk '$1!=""'  > "${genomica2%.nii.gz}.ic$((i+1))"
	rm "${genomica1%.nii.gz}.txt$(printf %05d $i)" "${genomica2%.nii.gz}.txt$(printf %05d $i)"
done


cat $catalogue|awk '$2==1{print $1}'|sed 's/_/ /g' > ${catalogue}.5million
