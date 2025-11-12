#module load 2025 ; module load R/4.4.2-gfbf-2025a

library(dplyr)
library(tidyr)
library(ggplot2)

df <-read.table('/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats/pheno_psy_gencor.rho',header=T)
##df <-read.table('~/genetic-data/genica/gwas-databases/panukbb/sumstats/pheno_psy_gencor.z',header=T)

df_long <- df %>%
  pivot_longer(
    cols = -1,               # first col is phenotype
    names_to = "biomarker",
    values_to = "value"
  )

ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", na.rm = TRUE) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~ biomarker, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

