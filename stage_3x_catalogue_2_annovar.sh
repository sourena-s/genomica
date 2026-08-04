FSLDIR=/projects/0/einf2700/sourena/software/fsl
export FSLDIR
. ${FSLDIR}/etc/fslconf/fsl.sh

annovar_path=/projects/0/einf2700/sourena/software/annovar

#The below catalogue is the same for r2psy and manual (it's a merged mask)
catalogue=/home/ssoheili/genomica-code/clumps/merged_manual-r2point3_psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue

cat $catalogue|awk '$2==1{print $1}'|sed 's/_/ /g'|awk '{print $1" "$2" "$2+length($3)-1" "$3" "$4" "$1"_"$2"_"$3"_"$4}' > ${catalogue}.annovar   # $4 is the reference allele

catalogue=${catalogue}.annovar
genomic_ics=$1
genomic_ics=${genomic_ics%.nii.gz}

echo reading file $catalogue and $genomic_ics

ic_count=$($FSLDIR/bin/fslinfo $genomic_ics |awk '$1=="dim4"{print $2}')
var_count=$($FSLDIR/bin/fslinfo $genomic_ics |awk '$1=="dim1"{print $2}')

if [ ! $var_count -eq $(cat $catalogue|wc -l) ];then echo "Catalogue and ICA component variant count mismatch: $var_count vs $(cat $catalogue|wc -l)";exit 1;fi

$FSLDIR/bin/fsl2ascii $genomic_ics ${genomic_ics}_txt

for ((i=0;i<ic_count;i++));do
        ic_file_name=${genomic_ics}_txt$(printf %05d $i)

cat $ic_file_name |tr " " "\n"|awk '$1!=""' > ${ic_file_name}_rows
rm $ic_file_name
ic_file_name=${ic_file_name}_rows


echo "pasting to ${genomic_ics}_$((i+1))"
paste -d " " ${catalogue} $ic_file_name  > "${genomic_ics}_$((i+1))"


#cat "$stat_file" |awk '$8>7.3'|sort -k8,8 -n -r |awk '{print $2" "$3" "$3+length($5)-1" "$5" "$6" "$4" "$1" "$8}' > sig_$(echo $i|sed 's/\//_/g').annovar

$annovar_path/annotate_variation.pl -geneanno -buildver hg19 -dbtype refgene "${genomic_ics}_$((i+1))" $annovar_path/humandb/

cat ${genomic_ics}_$((i+1)).variant_function|sed 's/\t/ /g' |gzip > ${genomic_ics}_$((i+1)).variant_function_ssv.tar.gz

rm ${genomic_ics}_$((i+1)) $ic_file_name

done


