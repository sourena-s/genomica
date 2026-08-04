#module load 2024; module load R/4.4.2-gfbf-2024a

library('bigsnpr')
library(data.table)
base_path <- '/projects/0/einf2700/sourena/genica-2025/gwas-databases'
snp_table <- paste0(base_path, '/panukbb/full_variant_qc_metrics.sorted.sst.rsid')
psy_list <- read.table('psy_gwas_list.txt')[,1]

info_snp <-fread(snp_table,header=T)

for (gwas_file in psy_list) {

print (paste('Reading GWAS sumstats:', gwas_file))
output_file <- paste0(gwas_file, '.ukb.matched')

gwas <- fread(gwas_file,header=T)
gwas$chr <- as.character(gwas$chr)
#change z column name to beta

#remove beta and keep z for matching
gwas <- subset(gwas,select = c("rsid" ,"chr", "pos", "a1", "a0", "p", "z"))
gwas$beta <- gwas$z
matched_gwas <- snp_match(sumstats=gwas,info_snp=info_snp, strand_flip=FALSE, return_flip_and_rev = TRUE)
matched_gwas$matched_z <- matched_gwas$beta

matched_gwas$matched_a1 <- ifelse(matched_gwas[["_REV_"]], matched_gwas$a0, matched_gwas$a1)
matched_gwas$matched_a0 <- ifelse(matched_gwas[["_REV_"]], matched_gwas$a1, matched_gwas$a0)

matched_gwas <- subset(matched_gwas,select = c("chr", "pos", "a0", "a1", "rsid" ,"rsid.ss" , "_REV_" , "matched_a0", "matched_a1", "p", "z", "matched_z"))
write.table(matched_gwas, file = output_file, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
