module load 2025; module load Anaconda3/2025.06-1;source activate ldsc

ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py
sumstat=/home/ssoheili/genetic-data/genica/gwas-databases/PSY/ADHD/DEMONTIS2023/ADHD2022_iPSYCH_deCODE_PGC.meta.with_logodds.matched
ref_path=/home/ssoheili/genetic-data/genica/software/ldsc/ref

header_file=${sumstat%*matched}header

n_tot=$(awk '$1=="n_tot"{print $2}' $header_file)
if [ $n_tot == ""];then n_str="--N-cas-col ncas --N-con-col ncon"; else n_str="--N $n_tot";fi

$munge_bin --sumstats $sumstat $n_str --info-min 0.8 --snp rsid --a1 a0 --a2 a1 --p p --signed-sumstat beta_matched,0 --info impinfo --a1-inc --out ${sumstat}.ldsc ### --maf-min 0.01 --frq

##$ldsc_bin --h2 ${sumstat}.ldsc.sumstats.gz --ref-ld $ref_path/ukbb/UKBB.EUR.rsid \
##                --out "${sumstat}.ldsc.UKBB.h2" \
##		--frqfile-chr $ref_path/1000GP$ref_path/1000GP/1000G_Phase3_frq/1000G.EUR.QC. \
##                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

##$ldsc_bin --h2 ${sumstat}.ldsc.sumstats.gz --ref-ld $ref_path/ukbb/UKBB.EUR.25LDMS.rsid \
##                --out "${sumstat}.ldsc.UKBB.h2" \
##		--frqfile-chr $ref_path/1000GP$ref_path/1000GP/1000G_Phase3_frq/1000G.EUR.QC. \
##                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

##$ldsc_bin --h2 ${sumstat}.ldsc.sumstats.gz --ref-ld $ref_path/ukbb/UKBB.EUR.8LDMS.rsid \
##                --out "${sumstat}.ldsc.UKBB.h2" \
##		--frqfile-chr $ref_path/1000GP$ref_path/1000GP/1000G_Phase3_frq/1000G.EUR.QC. \
##                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

$ldsc_bin --h2 ${sumstat}.ldsc.sumstats.gz --ref-ld-chr $ref_path/1000GP/baselineLD. \
                --out "${sumstat}.ldsc.1000GP.h2" \
                --w-ld-chr $ref_path/1000GP/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

