#!/usr/bin/env bash
input_gwas_list="psy_gwas_list.txt"

module load 2024
module load parallel/20240722-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0

python_bin=python3 

var_catalogue=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted
var_catalogue_rsid=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst.rsid.txt

#dl_path=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats
dl_path=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats-psy

cat $input_gwas_list|while read mat_file;do

	mat_file="${mat_file}.ukb.matched"
	echo "processing gwas $mat_file"

    header_file="${mat_file}.header"
    
    head -n 1 "${mat_file}" |tr " " "\n" > "$header_file"
    echo "nrows $(cat "${mat_file}" |wc -l)" >> "$header_file"
    awk -F' ' 'NR>1{print $1"_"$2"_"$3"_"$4" "$0}' "${mat_file}" |LC_ALL=C sort -t' ' -k1,1 > "${mat_file}.sorted"
    LC_ALL=C join -a1 -j1 "$var_catalogue" "${mat_file}.sorted" > "${mat_file}.ref" #in reference catalogue table
    cat "${mat_file}.sorted" |cut -d' ' -f1 > "${mat_file}.sorted.vars"

z_col=$(awk -v name="matched_z" '$1==name {print NR}' "$header_file")
p_col=$(awk -v name="p" '$1==name {print NR}' "$header_file")

        echo "Processing z_col=$z_col +1, p_col=$p_col"

	awk -F' ' -v z_col=$((z_col+1)) -v p_col=$((p_col+1)) '
        { cond=($z_col!="NA" && $p_col!="NA" && $p_col>=0)
                if (cond) {
                    z_stat = $z_col
                } else {z_stat = 0 }
	print z_stat
            }
        ' "${mat_file}.ref" > "${mat_file}.z"

	####run ldsc munge and calculate corr with PSY GWASes
	echo "Processing ldsc: p_col=${p_col} +1"
	awk -F' ' -v p=$((p_col+1)) '{print $p}' "${mat_file}.ref" > "${mat_file}.p"

	echo "Munging ${mat_file}.z and variant file ${mat_file}.sorted.vars"
	echo "chr pos a0 a1 rsid z p" > "${mat_file}.tomunge"
	paste -d" " "$var_catalogue_rsid" "${mat_file}.z" "${mat_file}.p" |awk '$6!="nan" && $7!="NA"' >> "${mat_file}.tomunge"
#	module load 2025; module load Anaconda3/2025.06-1;source activate ldsc
	ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
	munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py
	ref_path=/projects/0/einf2700/sourena/software/ldsc/ref
	psy_path=/projects/0/einf2700/sourena/genica-2025/gwas-databases/PSY
	py2_path=/home/ssoheili/.conda/envs/ldsc/bin/python
	n_str="--N 10000"
	
	$py2_path	$munge_bin --sumstats "${mat_file}.tomunge" $n_str --chunksize 500000 --merge-alleles ${ref_path}/w_hm3.snplist --p p --snp rsid --a1 a1 --a2 a0 --signed-sumstat z,0 --out "${mat_file}.ldsc"

	echo "Checking if munged sumstat file exists:"
	ls ${mat_file}.ldsc.sumstats.gz
	echo "Proceeding to rG estimation.."
	psy_array=("AD" "ASD" "ANX" "BIP_meta_2025" "BIP_clinical" "BIP_community" "MDD" "ADHD" "SCZ" "OCD" "PTSD" "ADDICTION")
	psy_gwases=$(for psy in "${psy_array[@]}";do echo -n ",${psy_path}/$psy/gwas.orig.sumstats.gz";done)
	
	$py2_path $ldsc_bin --h2 "${mat_file}.ldsc.sumstats.gz" --ref-ld-chr $ref_path/1000GP/global/LDscore/LDscore. \
	--frqfile-chr $ref_path/1000GP/global/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --out "${mat_file}.ldsc.h2"

	$py2_path $ldsc_bin --rg "${mat_file}.ldsc.sumstats.gz$psy_gwases" --ref-ld-chr $ref_path/1000GP/global/LDscore/LDscore. \
	--frqfile-chr $ref_path/1000GP/global/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --out "${mat_file}.ldsc.psy.corr"

        awk -vs=${mat_file}.ldsc.sumstats.gz '$1==s' ${mat_file}.ldsc.psy.corr.log | sed 's/.EUR.ldsc.sumstats.gz//g'| awk '{$1=$1}1' OFS=' ' > ${mat_file}.psy_rg

	echo "Phenotype ${psy_array[@]}" > ${mat_file}.psy_rg.rho
	cp ${mat_file}.psy_rg.rho ${mat_file}.psy_rg.z
	echo -e "${mat_file} $(awk -F' ' '{print  $3}' "${mat_file}.psy_rg"|tr "\n" " ")\n" >> ${mat_file}.psy_rg.rho
	echo -e "${mat_file} $(awk -F' ' '{print $5}' "${mat_file}.psy_rg"|tr "\n" " ")\n" >> ${mat_file}.psy_rg.z

#    $python_bin /home/ssoheili/genomica-code/panukb_dl_np_compress.py "${mat_file}.z"

#     rm -rf ${mat_file}.tomunge ${mat_file}.p ${mat_file}.ldsc.sumstats.gz ${mat_file}.sorted.vars ${mat_file}.ref ${mat_file}.sorted ## ${mat_file}.z 
 
done



clump_mask=/home/ssoheili/genomica-code/clumps/merged_manual-r2point3_psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue
clump_mask_size=$(awk '$2==1' $clump_mask|wc -l)
export FSLOUTPUTTYPE="NIFTI2_GZ"
ica1=/projects/0/einf2700/sourena/genica-2025/melodic_results/matrix_4d_manual_mergedmask_mergedmask_ICA_dim_20_nl_pow3/melodic_oIC.nii.gz
ica2=/projects/0/einf2700/sourena/genica-2025/melodic_results/matrix_4d_r2psy_mergedmask_mergedmask_ICA_dim_20_nl_pow3/melodic_oIC.nii.gz

cat $input_gwas_list|while read s;do
echo processing $s
paste -d" " $clump_mask ${s}.ukb.matched.z |awk '$2==1{print $5}' > ${s}.ukb.matched.z.clumped
paste -d" " $clump_mask ${s}.ukb.matched.ref|awk '$2==1{if ($6=="") {print 0;}else{print 1;}}' > ${s}.ukb.matched.z.clump.overlap_mask

/projects/0/einf2700/sourena/software/fsl/bin/fslascii2img ${s}.ukb.matched.z.clumped $clump_mask_size 1 1 1 1 1 1 1 ${s}.ukb.matched.z.nifti
/projects/0/einf2700/sourena/software/fsl/bin/fslascii2img ${s}.ukb.matched.z.clump.overlap_mask $clump_mask_size 1 1 1 1 1 1 1 ${s}.ukb.matched.z.clump_overlap_mask.nifti

cat ${s}.ukb.matched.z.clump.overlap_mask|grep 1|wc -l > ${s}.ukb.matched.z.clump.overlap_mask.size

/projects/0/einf2700/sourena/software/fsl/bin/fsl_glm -i ${s}.ukb.matched.z.nifti -m ${s}.ukb.matched.z.clump_overlap_mask.nifti -d $ica1 -o ${s}.ukb.matched.z.nifti.ica1_mix
/projects/0/einf2700/sourena/software/fsl/bin/fsl_glm -i ${s}.ukb.matched.z.nifti -m ${s}.ukb.matched.z.clump_overlap_mask.nifti -d $ica2 -o ${s}.ukb.matched.z.nifti.ica2_mix

done

echo "condition overlap_mask_nvar $(for i in manual psy_r2;do for ((j=1;j<=20;j++));do echo "${i}_ic_$j";done;done|tr "\n" " ")" >  psy_gwas_mix

cat psy_gwas_list.txt |awk -F'/' '{print $9}' > psy_gwas_name
cat $input_gwas_list|while read s;do paste -d" " ${s}.ukb.matched.z.clump.overlap_mask.size ${s}.ukb.matched.z.nifti.ica1_mix ${s}.ukb.matched.z.nifti.ica2_mix;done > paste_psy_ica1_ica2

paste -d" " psy_gwas_name paste_psy_ica1_ica2 >> psy_gwas_mix

rm psy_gwas_name paste_psy_ica1_ica2
