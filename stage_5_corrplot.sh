# module load 2025 ; module load R/4.4.2-gfbf-2025a

dims=(5 10 15 20)

#for dim1 in "${dims[@]}";do for dim2 in "${dims[@]}";do  fslcc -p 6 --noabs -t -2 --nodemean matrix_4d_manual_rpoint3_ICA_dim_$dim1/melodic_oIC.nii.gz matrix_4d_manual_rpoint3_ICA_dim_$dim2/melodic_oIC.nii.gz;done;done > corr_man2man.txt
#for dim1 in "${dims[@]}";do for dim2 in "${dims[@]}";do  fslcc -p 6 --noabs -t -2 --nodemean matrix_4d_r2psy_ICA_dim_$dim1/melodic_oIC.nii.gz matrix_4d_r2psy_ICA_dim_$dim2/melodic_oIC.nii.gz;done;done > corr_psy2psy.txt



for melodic_str in matrix_4d_manual_mergedmask_mergedmask_ICA_dim_ matrix_4d_r2psy_mergedmask_mergedmask_ICA_dim_;do
	for nl in _nl_pow3  _nl_pow4  _nl_tanh;do

for dim1 in "${dims[@]}";do for dim2 in "${dims[@]}";do  fslcc -p 6 --noabs -t -2 --nodemean matrix_4d_r2psy_mergedmask_mergedmask_ICA_dim_${dim1}_nl_pow3/melodic_oIC.nii.gz matrix_4d_r2psy_mergedmask_mergedmask_ICA_dim_${dim2}_nl_pow3/melodic_oIC.nii.gz;done;done > corr_man2man_${melodic_str}_${nl_str}.txt
for dim1 in "${dims[@]}";do for dim2 in "${dims[@]}";do  fslcc -p 6 --noabs -t -2 --nodemean matrix_4d_manual_mergedmask_mergedmask_ICA_dim_${dim1}_nl_pow3/melodic_oIC.nii.gz matrix_4d_manual_mergedmask_mergedmask_ICA_dim_${dim2}_nl_pow3/melodic_oIC.nii.gz;done;done > corr_psy2psy.txt

ic1_loc=0
ic2_loc=0

for dim1 in "${dims[@]}"
	do for ((ic1=1;ic1<=dim1;ic1++))
	do
	ic1_loc=$((ic1_loc+1))
	ic2_loc=0
		for dim2 in "${dims[@]}"
			do for ((ic2=1;ic2<=dim2;ic2++))
			do
				ic2_loc=$((ic2_loc+1))
				echo "Dim-${dim1}-ic-${ic1}-Dim-${dim2}-ic-${ic2} $ic1_loc $ic2_loc"
			done
		done
	done
done |sort -k1,1 > ic_labels.ids




for dim1 in "${dims[@]}"
	do for dim2 in "${dims[@]}"
		do for ((ic1=1;ic1<=dim1;ic1++))
			do for ((ic2=1;ic2<=dim2;ic2++))
			do echo "Dim-${dim1}-ic-${ic1}-Dim-${dim2}-ic-${ic2}"

done;done;done;done > ic_labels

paste -d" " ic_labels corr_psy2psy.txt |sort -k1,1 > corr_psy2psy.1

join -j1 ic_labels.ids corr_psy2psy.1 |	sort -k2,2n -k3,3n |cut -d" " -f1-3,6 >  corr_psy2psy.2


paste -d" " ic_labels corr_man2man.txt |sort -k1,1 > corr_man2man.1

join -j1 ic_labels.ids corr_man2man.1 |	sort -k2,2n -k3,3n |cut -d" " -f1-3,6 >  corr_man2man.2


#for dim1 in "${dims[@]}";do for dim2 in "${dims[@]}";do  fslcc -p 6 --noabs -t -2 --nodemean matrix_4d_r2psy_ICA_dim_$dim1/melodic_oIC.nii.gz matrix_4d_manual_rpoint3_ICA_dim_$dim2/melodic_oIC.nii.gz;done;done > corr_psy2man.txt


