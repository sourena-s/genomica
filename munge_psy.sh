#!/usr/bin/env bash
module load 2025; module load Anaconda3/2025.06-1;source activate ldsc

ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py
base_dir=/home/ssoheili/genetic-data/genica/gwas-databases/PSY
ref_path=/home/ssoheili/genetic-data/genica/software/ldsc/ref

info_str="--info impinfo" 

for trait in PTSD ASD ADDICTION ANX BIP OCD MDD ADHD SCZ;do

case "$trait" in
	"ADDICTION") n_str="--N 1025550"
		;;
	"ANX") n_str="--N 1096458"
		;;
	"ASD")  n_str="--N-cas 18382 --N-con 27969"
		;;
	"PTSD") n_str="--N-col n"
		info_str=""
		;;

	*) n_str="--N-cas-col ncas --N-con-col ncon"
		;;
esac

	
##$munge_bin --sumstats "$base_dir/$trait/gwas.matched" $n_str $info_str --p p --snp rsid --signed-sumstat matched_beta,0 --out "$base_dir/$trait/gwas.matched.ldsc" --merge-alleles ${ref_path}/w_hm3.snplist  --a1 a1 --a2 a0 --chunksize 10000

case "$trait" in
	"PTSD") signed_stat="z" 
		;;
	*) signed_stat="beta" 
		;;
esac
	

# $munge_bin --sumstats "$base_dir/$trait/gwas.orig" $n_str $info_str --p p --snp rsid --signed-sumstat ${signed_stat},0 --out "$base_dir/$trait/gwas.orig" --merge-alleles ${ref_path}/w_hm3.snplist  --a1 a1 --a2 a0 --chunksize 10000

###### $ldsc_bin --h2 $base_dir/$trait/gwas.matched.ldsc.sumstats.gz --ref-ld-chr $ref_path/1000GP/baselineLD. --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --out $base_dir/$trait/gwas.matched.h2.1000GP
## $ldsc_bin --h2 $base_dir/$trait/gwas.matched.ldsc.sumstats.gz --ref-ld $ref_path/ukbb/UKBB.EUR.rsid --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --out $base_dir/$trait/gwas.matched.h2.UKB

###### $ldsc_bin --h2 ~/genetic-data/genica/gwas-databases/PSY/BIP/gwas.matched.ldsc.sumstats.gz --ref-ld $ref_path/ukbb/UKBB.EUR.rsid --out test --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --out test
###$ldsc_bin --h2 "$base_dir/$trait/gwas.matched.ldsc" --ref-ld $ref_path/UKBB.EUR.rsid --out $base_dir/$trait/gwas.matched.h2
	##--w-ld $ref_path/eur_w_ld_chr/ \
	##--overlap-annot UKBB.EUR.25LDMS \

done



 psy_path=/home/ssoheili/genetic-data/genica/gwas-databases/PSY
psy_array=("ADHD" "SCZ" "ASD" "ANX" "BIP" "MDD" "OCD" "PTSD" "ADDICTION")
        psy_gwases=$(for psy in "${psy_array[@]}";do echo -n "${psy_path}/$psy/gwas.orig.sumstats.gz,";done)

        $ldsc_bin --rg "${psy_gwases%*,}" --ref-ld-chr $ref_path/1000GP/baselineLD. \
                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                --out "estimated_psy_corr"

