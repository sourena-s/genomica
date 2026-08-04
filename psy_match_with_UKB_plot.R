
librar(tidyverse)
library(ggplot)

a<-read.table('psy_gwas_mix', header=T)

a <- a %>% mutate(GWAS_label = paste0(condition, " (nvar=", overlap_mask_nvar, ")"))

component_order <- setdiff(colnames(a),c("GWAS_label","condition","overlap_mask_nvar"))
row_order <- a$GWAS_label

long_df <- a %>% pivot_longer(cols = all_of(component_order), names_to = "Component", values_to = "Value") %>%
    mutate(Component = factor(Component, levels = component_order),
    GWAS_label = factor(GWAS_label, levels = unique(GWAS_label))) %>%
mutate(GWAS_label = factor(GWAS_label, levels = row_order))

max_abs <- max(abs(long_df$Value))

ggplot(long_df, aes(Component, GWAS_label, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-max_abs, max_abs)
  ) + theme_bw() + theme(axis.text.x = element_text(angle = 90))

