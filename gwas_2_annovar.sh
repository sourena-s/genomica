#assumes the first allele ($4) to be the reference. Need to check with UKB genotype format

#ID CHROM POS ORIG_RSID REF ALT A1 NEG_LOG10_P
#1:756604[call]A,G 1 756604 rs3131962 A G A 2.98305
#1:761147[call]T,C 1 761147 rs3115850 T C T 2.87705
#1:767096[call]A,G 1 767096 rs115991721 A G G 2.99459
#1:768448[call]G,A 1 768448 rs12562034 G A A 2.42034
#1:779322[call]A,G 1 779322 rs4040617 A G G 2.87987


10_100000012_G_A 0
10_100000056_C_G 1
10_10000011_G_A 1
10_100000122_T_A 0
10_10000018_A_G 0
10_10000030_A_C 0
10_10000033_G_A 0
10_100000588_T_C 0



catalogue=$1
genomic_ics=$2
output_file=$3

echo reading file $catalogue and $genomic_ics

awk '$2==1{print $1}' $catalogue > ${catalogue}_vars_in_mask

ic_count=$(fslinfo $genomic_ics |awk '$1=="dim4"{print $2}')
var_count=$(fslinfo $genomic_ics |awk '$1=="dim1"{print $2}')

if [ ! $var_count -eq $(awk '$2==1' $catalogue|wc -l) ];then echo "Catalogue and ICA component variant count mismatch: $var_count vs $(awk '$2==1' $catalogue|wc -l)";exit 1;fi

fsl2ascii $genomic_ics ${genomic_ics}_txt

for ((i=0;i<ic_count;i++));do
	ic_file_name=${genomic_ics}_txt$(printf %05d $i)

cat $ic_file_name |tr " " "\n"|awk '$1!=""' > ${ic_file_name}_rows
ic_file_name=${ic_file_name}_rows
paste -d " " ${catalogue}_vars_in_mask $ic_file_name > "${genomic_ics}_$((i+1))"

cat "${genomic_ics}_$((i+1))" | sed 's/_/ /g'|awk '$1" "$2" "$2+length($3)-1" "$3" "$4" "rsid" "$1" "pval' > "${genomic_ics}_$((i+1))_annovar"

annovar=/home/ssoheili/genetic-data/genica/software/annovar/annotate_variation.pl


#cat "$stat_file" |awk '$8>7.3'|sort -k8,8 -n -r |awk '{print $2" "$3" "$3+length($5)-1" "$5" "$6" "$4" "$1" "$8}' > sig_$(echo $i|sed 's/\//_/g').annovar

$annovar -geneanno -buildver hg19 -dbtype refgene $genomic_ics_$((i+1))_annovar ~/software/annovar/humandb/

done


