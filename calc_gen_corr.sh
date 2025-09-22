module load 2025; module load Anaconda3/2025.06-1;source activate ldsc
ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py
sumstat=/home/ssoheili/genetic-data/genica/gwas-databases/PSY/ADHD/DEMONTIS2023/ADHD2022_iPSYCH_deCODE_PGC.meta.with_logodds.matched
$munge_bin --sumstats $sumstat --N-cas-col nca --N-con-col nco --info-min 0.3 --snp rsid --a1 a0 --a2 a1 --p p --signed-sumstat beta_matched,0 --info impinfo --a1-inc --out ${sumstat}.ldsc ### --maf-min 0.01 --frq --N 
