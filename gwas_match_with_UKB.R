##module load 2024; module load R/4.4.2-gfbf-2024a

library('bigsnpr')
#snp_table <- '/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst'
base_path <- '/projects/0/einf2700/sourena/genica-2025/gwas-databases'

snp_table <- paste0(base_path, '/panukbb/full_variant_qc_metrics.sorted.sst.rsid')

gwas_list <- c(
	  	'ADDICTION/Hatoum2023AddictionEuropean', 
	        'ASD/iPSYCH-PGC_ASD_Nov2017_with_logodds',
		'SCZ/Trubetskoy2022/PGC3_SCZ_wave3.european.autosome.public.v3',
		'ADHD/DEMONTIS2023/ADHD2022_iPSYCH_deCODE_PGC.meta.with_logodds',
		'BIP_meta_2025/Oconnell2025/bip2024_eur_no23andMe.with_logodds',
		'BIP_clinical/gwas',
		'BIP_community/gwas',
		'ANX/Lokhammer2024/ANX_EUR',
		'MDD/CELL2025/pgc-mdd2025_Clin_eur_v3-49-24-11',
		'OCD/STROM2024/ocs2024obsessive-compulsive_symptoms_daner_STR_NTR_SfS_TwinsUK_strometal.with_logodds',
	  	'PTSD/Nievergelt2024/eur_ptsd_pcs_v4_aug3_2021',
		'AD/Wightman_2021/PGCALZ2sumstatsExcluding23andMe.rsid'
)

gwas_list <- paste0(base_path,'/PSY/',gwas_list)

for (gwas_filepath in gwas_list) {

print (paste('Reading GWAS sumstats:', gwas_filepath))
gwas <- paste0(gwas_filepath,'.txt')
gwas_header <- paste0(gwas_filepath,'.header')
output_file <- paste0(gwas_filepath, '.matched')
output_file_orig <- paste0(gwas_filepath, '.orig')

info_snp <-read.table(snp_table,header=T)

gwas <-read.table(gwas,header=T)
gwas_header <-read.table(gwas_header,header=F)[,2]

colnames(gwas) <- gwas_header[1:ncol(gwas)] #last rows in .header file may be n_tot or other info
gwas$chr <- as.character(gwas$chr)

#matched_gwas <- snp_match(sumstats=gwas,info_snp=info_snp, strand_flip=FALSE, return_flip_and_rev = TRUE)

## matched_gwas$matched_OR <- exp(matched_gwas$beta)

## if ("se" %in% colnames(matched_gwas)) {matched_gwas$matched_z <- matched_gwas$beta / matched_gwas$se }
## if ("or" %in% colnames(matched_gwas)) {names(matched_gwas)[names(matched_gwas) == "or"] <- "orig_OR" }

## names(matched_gwas)[names(matched_gwas) == "beta"] <- "matched_beta"

#write.table(matched_gwas, file = output_file, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(gwas, file = output_file_orig, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
