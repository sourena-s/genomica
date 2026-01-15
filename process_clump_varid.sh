#converts clump variant IDs into strings for joining with h5 variant catalogue
cat clumps.all |sed 's/\[[^]]*\]/_/g'|sed 's/:/_/g'|sed 's/,/_/g'|LC_ALL=C sort -k1,1 |awk '{print $1" 1"}'  > clumps.sorted

LC_ALL=C join -a1 -j1 /projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted clumps.sorted |awk '{if ($2=="") $2=0;print $0}' > clumps.sorted.mask.in.var.catalogue
