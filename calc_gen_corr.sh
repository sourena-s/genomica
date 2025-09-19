module load 2025; module load Anaconda3/2025.06-1
ldsc_bin=/home/ssoheili/software/ldsc/ldsc.py
munge_bin=/home/ssoheili/software/ldsc/munge_sumstats.py

$munge_bin --sumstats --N --N-cas --N-con --info-min --maf-min --snp --a1 a0 --a2 a1 --p p --frq --signed-sumstat beta --info impinfo --a1-inc
