#!/usr/bin/env bash
input_gwas_list=/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/gwas_download_table.reformat.B
input_gwas_list=/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/phenotype_manifest.proc.B 
input_gwas_list=/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/phenotype_manifest.proc.B.16test 

NPROC=$2
if [ $NPROC == "" ];then NPROC=`nproc --all`;fi
echo "Running on $NPROC cores"

array_id=$1

start_row=$(echo "( $array_id * $NPROC ) + 1"|bc) #inclusief
end_row=$(echo "( $array_id + 1) * $NPROC"|bc) #inclusief

echo "Received job array index $array_id corresponding to rows >= $start_row and < $end_row"

export OMP_PLACES cores
export OMP_PROC_BIND true

module load 2024
module load parallel/20240722-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0
## module load 2025; module load Anaconda3/2025.06-1

# download, unzip, process
process_url() {
python_bin=python3 # which python? use anaconda?
dl_path=/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/sumstats
var_catalogue=/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/full_variant_qc_metrics.sorted
var_catalogue_rsid=/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst.rsid.txt
if [ ! -d "$dl_path/gwas_meta" ];then mkdir -p "$dl_path/gwas_meta";fi

    url=$1
    ncases=$2
    ncontrols=$3
    mat_file=${url##*sumstats_flat_files/}
    mat_file=${mat_file%.tsv.bgz*}

    echo "Source URL: $url" 
    echo "Saving to ${mat_file}.dl"
    wget -q "$url" -O "$dl_path/gwas_meta/${mat_file}.dl"

    header_file="$dl_path/gwas_meta/${mat_file}.header"
    
    zcat "$dl_path/gwas_meta/${mat_file}.dl" > "$dl_path/gwas_meta/${mat_file}.orig" #unsorted
    head -n 1 "$dl_path/gwas_meta/${mat_file}.orig" |tr "\t" "\n" > "$header_file"
    echo "nrows $(cat "$dl_path/gwas_meta/${mat_file}.orig" |wc -l)" >> "$header_file"
    awk -F$'\t' 'NR>1{print $1"_"$2"_"$3"_"$4"\t"$0}' "$dl_path/gwas_meta/${mat_file}.orig" |LC_ALL=C sort -t$'\t' -k1,1 > "$dl_path/gwas_meta/${mat_file}.sorted"
    join -a1 -j1 "$var_catalogue" "$dl_path/gwas_meta/${mat_file}.sorted" > "$dl_path/gwas_meta/${mat_file}.ref" #in reference catalogue table
    cat "$dl_path/gwas_meta/${mat_file}.sorted" |cut -d$'\t' -f1 > "$dl_path/gwas_meta/${mat_file}.sorted.vars"

populations=("meta" "meta_hq" "AFR" "CSA" "AMR" "EAS" "EUR" "MID")
populations=("EUR")

for pop in "${populations[@]}"; do
    beta_name="beta_$pop"
    se_name="se_$pop"

    beta_col=$(awk -v name="$beta_name" '$1==name {print NR}' "$header_file")
    se_col=$(awk -v name="$se_name" '$1==name {print NR}' "$header_file")
    eur_p_col=$(awk -v name="neglog10_pval_EUR" '$1==name {print NR}' "$header_file")

    if [ "$pop" == "meta" ] || [ "$pop" == "meta_hq" ]; then
        lowconf_opt=""
    else
       var_lowconf_col=$(awk -v name="low_confidence_$pop" '$1==name {print NR}' "$header_file")
       lowconf_opt="-v lowconf=$((var_lowconf_col+1))"
    fi


    if [ -n "$beta_col" ] && [ -n "$se_col" ]; then

	if [ ! -d "$dl_path/$pop" ];then echo "making dir $dl_path/$pop"; mkdir -p "$dl_path/$pop";fi

        echo "Processing $pop: beta_col=$beta_col +1, se_col=$se_col +1 and var QC awk argument $lowconf_opt"
	awk -F' ' -v b=$((beta_col+1)) -v s=$((se_col+1)) $lowconf_opt '
	{ cond=($s>0 && $s!="NA" && $b!="NA")
	if (lowconf!="") cond=(cond && ($lowconf!="true"))
		if (cond) {
                    z = $b/$s
                } else {z = "nan" }
                print z
            }
        ' "$dl_path/gwas_meta/${mat_file}.ref" > "$dl_path/$pop/${mat_file}.z.$pop"

	if [[ "$pop" == "EUR" ]];then 
	####run ldsc munge and calculate corr with PSY
	echo "Processing $pop: eur_p_col=$eur_p_col +1"
	awk -F' ' -v p=$((eur_p_col+1)) '{print $p}' "$dl_path/gwas_meta/${mat_file}.ref" > "$dl_path/$pop/${mat_file}.p.$pop"

	echo "Munging $dl_path/$pop/${mat_file}.z.$pop and variant file $dl_path/gwas_meta/${mat_file}.sorted.vars"
	echo "chr pos a0 a1 rsid z p" > "$dl_path/$pop/${mat_file}.${pop}.tomunge"
	paste -d" " "$var_catalogue_rsid" "$dl_path/$pop/${mat_file}.z.$pop" "$dl_path/$pop/${mat_file}.p.$pop" |awk '$6!="nan" && $7!="NA"{$7=(10 ** -($7));print $0;}' >> "$dl_path/$pop/${mat_file}.${pop}.tomunge"
	module load 2025; module load Anaconda3/2025.06-1;source activate ldsc
	ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
	munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py
	ref_path=/home/ssoheili/genetic-data/genica/software/ldsc/ref
	psy_path=/home/ssoheili/genetic-data/genica/gwas-databases/PSY
	
	if [[ $ncontrols == "NA" ]];then n_str="--N $ncases";else n_str="--N-cas $ncases --N-con $ncontrols";fi
	
	$munge_bin --sumstats "$dl_path/$pop/${mat_file}.${pop}.tomunge" $n_str --chunksize 10000 --merge-alleles ${ref_path}/w_hm3.snplist --p p --snp rsid --a1 a1 --a2 a0 --signed-sumstat z,0 --a1-inc --out "$dl_path/$pop/${mat_file}.${pop}.ldsc"

	echo "Checking if munged sumstat file exists:"
	ls $dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz
	echo "Proceeding to rG estimation.."
	psy_array=("ASD" "ANX" "BIP" "MDD" "ADHD" "SCZ" "OCD" "PTSD" "ADDICTION" )
	psy_gwases=$(for psy in "${psy_array[@]}";do echo -n ",${psy_path}/$psy/gwas.matched.ldsc.sumstats.gz";done)
	
	$ldsc_bin --rg "$dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz$psy_gwases" --ref-ld-chr $ref_path/1000GP/baselineLD. \
                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                --out "$dl_path/$pop/${mat_file}.${pop}.ldsc.psy.corr"

        awk -vs=$dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz '$1==s' $dl_path/$pop/${mat_file}.${pop}.ldsc.psy.corr.log | sed 's/.EUR.ldsc.sumstats.gz//g' > $dl_path/$pop/${mat_file}.${pop}.psy_rg

	awk '{printf "%s ", $3} END {print ""}' $dl_path/$pop/${mat_file}.${pop}.psy_rg|tr "\n" " "  > $dl_path/$pop/${mat_file}.${pop}.psy_rg.rho
	echo "${psy_array[@]}" > $dl_path/$pop/${mat_file}.${pop}.psy_rg.z
	awk '{printf "%s ", $5} END {print ""}' $dl_path/$pop/${mat_file}.${pop}.psy_rg|tr "\n" " " >> $dl_path/$pop/${mat_file}.${pop}.psy_rg.z
	fi

	#Save z as h5
    $python_bin /home/ssoheili/genomica-code/panukb_dl_np_compress.py "$dl_path/$pop/${mat_file}.z.$pop"
    
  ##  rm -rf $dl_path/$pop/${mat_file}.z.$pop $dl_path/$pop/${mat_file}.${pop}.tomunge $dl_path/$pop/${mat_file}.p.$pop $dl_path/$pop/${mat_file}.${pop}.ldsc.sumstats.gz

    fi
done
##rm -rf  $dl_path/gwas_meta/${mat_file}.sorted \
##	$dl_path/gwas_meta/${mat_file}.sorted.vars \
##	$dl_path/gwas_meta/${mat_file}.dl \
##	$dl_path/gwas_meta/${mat_file}.ref \
##	$dl_path/gwas_meta/${mat_file}.orig 
}
export -f process_url

# Feed URL and filename pairs from stdin to parallel
parallel -j "$NPROC" --colsep ' ' process_url {} :::: <(awk -F$'\t' 'NR>1{print $11"\t"$12"\t"$13}' "$input_gwas_list" |awk -F$'\t' -v start_row=$start_row -v end_row=$end_row 'NR >= start_row && NR <= end_row {print $1"\t"$2"\t"$3}'|while read filename ncas ncon;do echo https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/${filename} $ncas $ncon;done)

