#!/usr/bin/env bash

module load 2024
module load parallel/20240722-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0
## module load 2025; module load Anaconda3/2025.06-1

python_bin=python3 # which python? use anaconda?
dl_path=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats-psy
var_catalogue=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted
var_catalogue_rsid=/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst.rsid.txt

mat_file_list=psy_gwas_list.txt

cat $mat_file_list|while read s;do mat_file=${s}.ukb.matched
    #rsid chr pos a1 a0 p z matched_z

    awk -F' ' 'NR>1{print $2"_"$3"_"$5"_"$4" "$0}'  "$dl_path/gwas_meta/${mat_file}" |LC_ALL=C sort -t' ' -k1,1 > "$dl_path/gwas_meta/${mat_file}.var_id.sorted"
    echo "varid rsid chr pos a1 a0 p z matched_z" > "$dl_path/gwas_meta/${mat_file}.var_id.sorted.ref"
    LC_ALL=C join -a1 -j1 "$var_catalogue" "$dl_path/gwas_meta/${mat_file}.var_id.sorted" > "$dl_path/gwas_meta/${mat_file}.var_id.sorted.ref" #sumstats in ref


	####run ldsc munge and calculate corr with PSY GWASes
	awk -F' ' -v p=$((eur_p_col+1)) '{print $p}' "$dl_path/gwas_meta/${mat_file}.ref" > "$dl_path/$pop/${mat_file}.p.$pop"

	echo "Munging $dl_path/$pop/${mat_file}.z.$pop and variant file $dl_path/gwas_meta/${mat_file}.sorted.vars"
	echo "chr pos a0 a1 rsid z p" > "$dl_path/$pop/${mat_file}.${pop}.tomunge"
	paste -d" " "$var_catalogue_rsid" "$dl_path/$pop/${mat_file}.z.$pop" "$dl_path/$pop/${mat_file}.p.$pop" |awk '$6!="nan" && $7!="NA"{$7=(10 ** -($7));print $0;}' >> "$dl_path/$pop/${mat_file}.${pop}.tomunge"
#	module load 2025; module load Anaconda3/2025.06-1;source activate ldsc
	ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
	munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py
	ref_path=/projects/0/einf2700/sourena/software/ldsc/ref
	psy_path=/projects/0/einf2700/sourena/genica-2025/gwas-databases/PSY
	py2_path=/home/ssoheili/.conda/envs/ldsc/bin/python
	if [[ $ncontrols == "NA" ]];then n_str="--N $ncases";else n_str="--N-cas $ncases --N-con $ncontrols";fi
	
$py2_path	$munge_bin --sumstats "$dl_path/$pop/${mat_file}.${pop}.tomunge" $n_str --chunksize 500000 --merge-alleles ${ref_path}/w_hm3.snplist --p p --snp rsid --a1 a1 --a2 a0 --signed-sumstat z,0 --out "$dl_path/$pop/${mat_file}.${pop}.ldsc"

	echo "Checking if munged sumstat file exists:"
	ls $dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz
	echo "Proceeding to rG estimation.."
	psy_array=("AD" "ASD" "ANX" "BIP_meta_2025" "BIP_clinical" "BIP_community" "MDD" "ADHD" "SCZ" "OCD" "PTSD" "ADDICTION")
	psy_gwases=$(for psy in "${psy_array[@]}";do echo -n ",${psy_path}/$psy/gwas.orig.sumstats.gz";done)
	
##	$ldsc_bin --rg "$dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz$psy_gwases" --ref-ld-chr $ref_path/1000GP/baselineLD. \
##                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
##                --out "$dl_path/$pop/${mat_file}.${pop}.ldsc.psy.corr"

$py2_path $ldsc_bin --h2 "$dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz" --ref-ld-chr $ref_path/1000GP/global/LDscore/LDscore. \
	--frqfile-chr $ref_path/1000GP/global/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --out "$dl_path/$pop/${mat_file}.${pop}.ldsc.h2"

$py2_path $ldsc_bin --rg "$dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz$psy_gwases" --ref-ld-chr $ref_path/1000GP/global/LDscore/LDscore. \
	--frqfile-chr $ref_path/1000GP/global/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --out "$dl_path/$pop/${mat_file}.${pop}.ldsc.psy.corr"

        awk -vs=$dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz '$1==s' $dl_path/$pop/${mat_file}.${pop}.ldsc.psy.corr.log | sed 's/.EUR.ldsc.sumstats.gz//g'| awk '{$1=$1}1' OFS=' ' > $dl_path/$pop/${mat_file}.${pop}.psy_rg

	echo "Phenotype ${psy_array[@]}" > $dl_path/$pop/${mat_file}.${pop}.psy_rg.rho
	cp $dl_path/$pop/${mat_file}.${pop}.psy_rg.rho $dl_path/$pop/${mat_file}.${pop}.psy_rg.z
	echo -e "${mat_file} $(awk -F' ' '{print  $3}' "$dl_path/$pop/${mat_file}.${pop}.psy_rg"|tr "\n" " ")\n" >> $dl_path/$pop/${mat_file}.${pop}.psy_rg.rho
	echo -e "${mat_file} $(awk -F' ' '{print $5}' "$dl_path/$pop/${mat_file}.${pop}.psy_rg"|tr "\n" " ")\n" >> $dl_path/$pop/${mat_file}.${pop}.psy_rg.z
	fi

	#Save z as h5
    $python_bin /home/ssoheili/genomica-code/panukb_dl_np_compress.py "$dl_path/$pop/${mat_file}.z.$pop"


     rm -rf $dl_path/$pop/${mat_file}.${pop}.tomunge $dl_path/$pop/${mat_file}.p.$pop $dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz $dl_path/$pop/${mat_file}.z.$pop 

    fi
done
 rm -rf  $dl_path/gwas_meta/${mat_file}.sorted \
	$dl_path/gwas_meta/${mat_file}.sorted.vars \
	$dl_path/gwas_meta/${mat_file}.dl \
	$dl_path/gwas_meta/${mat_file}.ref \
	$dl_path/gwas_meta/${mat_file}.orig 
}
export -f process_url

# Feed URL and filename pairs from stdin to parallel
parallel -j "$NPROC" --colsep ' ' process_url {} :::: <(awk -F$'\t' 'NR>1{print $11"\t"$12"\t"$13}' "$input_gwas_list" |awk -F$'\t' -v start_row=$start_row -v end_row=$end_row 'NR >= start_row && NR <= end_row {print $1"\t"$2"\t"$3}'|while read filename ncas ncon;do echo https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/${filename} $ncas $ncon;done)

