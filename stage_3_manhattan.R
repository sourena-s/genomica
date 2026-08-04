#module load 2025 ; module load R/4.4.2-gfbf-2025a
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggnewscale)
library(svglite)
library(patchwork)
library(scales)
#stat_file= commandArgs(trailingOnly = T)[1]
#output_png= commandArgs(trailingOnly = T)[2]

plots<-list()
plot_title <- "Genomic IC"

#repel_thresh <- 8# ggrepel P-value threshold

y_min <- 0 
y_inter <- -log10((5)*10^-8)
#exonic_pval_threshold <- -log10((5)*10^-8)
  
catalogue <- '/home/ssoheili/genomica-code/clumps/merged_manual-r2point3_psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue'
catalogue <- fread(catalogue,header=F, sep = " ")
catalogue <- catalogue [V2==1,]$V1

#stat_file=paste0('/data/clusterfs/lag/users/sousoh/impact_modes/ica/', contrast, '/merged/main-alldim-minp.fil')

#ID CHROM POS ORIG_RSID REF ALT A1 NEG_LOG10_P
#10:94513[call]CCTGG,C 10 94513 Affx-52349171 CCTGG C C 0.176526

stat_file='matrix_4d_r2psy_mergedmask_mergedmask_ICA_dim_5_nl_pow3/melodic_oIC_1'

for (ic in seq(11,11)){

stat_file=paste0('matrix_4d_manual_mergedmask_mergedmask_ICA_dim_20_nl_pow3/melodic_oIC_',ic)

#intronic        R3HCC1L 10 100000056 100000056 G C 10_100000056_C_G -0.0688999
#intergenic      LOC101928272(dist=662455),LINC02670(dist=100674)        10 10000011 10000011 A G 10_10000011_G_A 0.145826

annovar1 <- fread( paste0(stat_file, '.variant_function.ssv') , header = FALSE, sep = " ")
annovar_cols <- c("Function", "Gene", "CHR", "Start", "End", "Ref", "Alt", "SNP", "IC_w")

#annovar_cols <- c("Function", "Gene", "CHR", "Start", "End", "Ref", "Alt", "SNP", "VariantTag", "P_annovar")
#annovar1 <- fread( paste0('sig_', contrast_str, '.annovar.variant_function.txt') , header = FALSE, sep = " ")
colnames(annovar1) <- annovar_cols
#annovar1 <- annovar1[annovar1$Function != "exonic",]
#annovar2 <-  fread( paste0(stat_file, '.exonic_variant_function') , header = FALSE)
#colnames(annovar2) <- annovar_cols
#annovar <-rbind(annovar1, annovar2)

#intronic        CPED1   7 120927437 120927437 A AC rs112901683 7:120927437[imp]A,AC 57.4682
#ID CHROM POS ORIG_RSID REF ALT A1 NEG_LOG10_P
#1:756604[call]A,G 1 756604 rs3131962 A G A 2.98305

print(stat_file)

a<-as.data.frame(annovar1)
#a$IC_w <- (a$IC_w - mean(a$IC_w)) / sd(a$IC_w)

m_z <- (a$IC_w - mean(a$IC_w)) / sd(a$IC_w - mean(a$IC_w))


m_p <- 2 * (1 - pnorm(abs(m_z)))
n <- length(m_p)

a$IC_w <- abs(a$IC_w / sd(a$IC_w))

pthresh_sparse <- 3.5
exonic_pval_threshold <- 4
repel_thresh <- 4
y_max <- max(a$IC_w)+ 2
top_chr_var_plot_pthresh <- 4

b<-a[!is.na(a$IC_w) & a$IC_w >= y_min,]

#pthresh_sparse <- -log10((5)*10^-8)

nbins <- 20
n_per_bin <- 200
print(dim(b))
b <- bind_rows(b %>% filter(IC_w > pthresh_sparse), b %>% filter(IC_w <= pthresh_sparse) %>% mutate(bin = cut(IC_w, breaks = nbins, include.lowest = TRUE)) %>% group_by(CHR,bin) %>% filter(n() > 0) %>% slice_sample(n = n_per_bin, replace = TRUE) %>% ungroup()) %>% select(-bin) ## %>% arrange(IC_w)

print(dim(b))
#colnames(b) <- c("VariantTag","CHR","BP","SNP","A1","A2","A3","P")

colnames(b) <- c("Function","Gene","CHR","BP","END_POS","A1","A2","SNP","IC_w")

######DEBUG
#qc_1<-as.data.frame(fread(qc_file,header=F,sep=" "))[,c(2,5)]
#qc_2<-as.data.frame(fread(qc_vartag,header=F,sep=" "))[,1]
#qc<-cbind(qc_1,qc_2)
#colnames(qc) <- c("SNP_QC","QC","VariantTag")
#print('Joining GWAS and QC list')
#m_joined <- left_join(b,qc,by="VariantTag")


#print('Merging GWAS and annotations')

#m_joined<- merge(m_joined, annovar[, c("VariantTag", "Gene", "Function")], by="VariantTag", all.x = TRUE)

#m_joined$Gene_OLD <- m_joined$Gene
#m_joined$Gene <- sub("\\(.*", "", m_joined$Gene_OLD)

m_joined <- b

top_exonic_snps_per_gene <- m_joined %>%
  filter(!is.na(Gene) & Function == "exonic" & IC_w > exonic_pval_threshold ) %>%
  group_by(Gene) %>%
  slice_max(order_by = IC_w, n = 1, with_ties = FALSE) %>% select(SNP)

print(paste('Found exonic variants', top_exonic_snps_per_gene))
 
top_snps_per_gene <- m_joined[m_joined$IC_w > repel_thresh,] %>%
  filter(!is.na(Gene) & Function != "intergenic") %>%
  group_by(Gene) %>%
  slice_max(order_by = IC_w, n = 1, with_ties = FALSE) %>% select(SNP) %>% filter(!SNP %in% top_exonic_snps_per_gene)


######################## d_f genesis
d_f <- m_joined %>%
  mutate(label = ifelse(SNP %in% top_snps_per_gene$SNP | SNP %in% top_exonic_snps_per_gene$SNP ,  paste0("'", SNP, "'~italic('(", Gene, ")')"), NA))

d_f$CHR[d_f$CHR %in% c("X", "XY")] <- 23
d_f$CHR[d_f$CHR == "Y"] <- 24
d_f$CHR[d_f$CHR == "MT"] <- 25

d_f$BP<-as.numeric(d_f$BP)
d_f$CHR <- as.numeric(d_f$CHR)

# Calculate cumulative base pair positions
offset <- 0
gap <- 20000000
d_f$dist <- 0

for (i in sort(unique(d_f$CHR)) ) {
	print(paste('chr ', i))
  idx <- (d_f$CHR == i) & !is.na(d_f$CHR) 
start <- min(d_f$BP[idx], na.rm = TRUE)
  stop <- max(d_f$BP[idx], na.rm = TRUE) 
  d_f$dist[idx] <- offset + d_f$BP[idx] - start
  offset <- offset + (stop - start) + gap
}

# Get unique chromosomes
d_f$CHR <- factor(d_f$CHR, levels = sort(unique(d_f$CHR)))
unique_chromosomes <- unique(d_f$CHR)

# Define the breaks and labels
chr_breaks <- tapply(d_f$dist, d_f$CHR, median)
rename_mapping <- c(1:22, "X", "Y", "MT")

# Assign names based on unique chromosomes, matching only as many as needed
chr_labels <- rename_mapping[1:length(unique_chromosomes)]
chr_labels <- rename_mapping[1:25]

#Sort by gene_id and descending P
#d_f_sorted <- d_f[order(d_f$Gene, -d_f$P), ]
#Keep only the first occurrence of each gene
#top_gene_ids <- d_f_sorted$gene_id[!duplicated(d_f_sorted$Gene)]
#In original data, set gene_id to NA if it's not in the top list
#d_f$gene_id[!(d_f$gene_id %in% top_gene_ids)] <- NA
d_f$gene_id <- NA
d_f$chr_max_var <- FALSE

for (i in sort(unique(d_f$CHR))) {
  chr_subset <- which(d_f$CHR == i & !is.na(d_f$IC_w) & d_f$IC_w > top_chr_var_plot_pthresh )
  if (length(chr_subset) > 0) {
    top_idx <- chr_subset[which.max(d_f$IC_w[chr_subset])]
    d_f$chr_max_var[top_idx] <- TRUE
    d_f$gene_id[top_idx] <- ifelse(d_f$Function[top_idx] == "intergenic", d_f$SNP[top_idx], paste0(d_f$SNP[top_idx], " ~italic('(", d_f$Gene[top_idx], ")')")  )	
#    d_f$gene_id[top_idx] <- ifelse(d_f$Function[top_idx] == "intergenic", NA, paste0(d_f$SNP[top_idx], " ~italic('(", d_f$Gene[top_idx], ")')")  )	
print(paste('Top CHR-', i, ' variants:' ,  d_f$gene_id[top_idx]))
}}

d_f <- d_f %>% mutate(gene_id = ifelse(SNP %in% top_exonic_snps_per_gene$SNP ,  label, gene_id) )
 
##one last round remove gene_id duplicates (e.g. imputed and called variants)
#d_f_sorted <- d_f[order(d_f$gene_id, -d_f$P), ]
#top_rows <- !duplicated(d_f_sorted$gene_id)
#d_f$gene_id[!(rownames(d_f) %in% rownames(d_f_sorted)[top_rows])] <- NA
# Calculate the length of the longest gene_id

max_length <- max(nchar(d_f$gene_id), na.rm = TRUE)



d_f$CHR <- as.numeric(d_f$CHR)

dpi=300 ##use 72 for SVG

d_f$gene_id <- gsub("_", ":", d_f$gene_id)

d_f <- subset(d_f, IC_w>0 & is.finite(dist) & is.finite(IC_w) ) 

d_f$var_type <- ifelse(d_f$Function %in% c("exonic"), "exonic", "nonexonic")

df_nonexonic <- subset(d_f,  ( Function != "exonic" ) | IC_w <= exonic_pval_threshold) 
df_exonic    <- subset(d_f,  IC_w > exonic_pval_threshold & !is.na(Function) & (Function == "exonic") |
		                 (chr_max_var &  (Function == "exonic")  )
		       ) 

chr_breaks <- chr_breaks[seq(1,25,1)]
chr_labels <- chr_labels[seq(1,25,1)]

d_f_unique_gene_id <- d_f %>%
  filter( !is.na(gene_id) &
 ( chr_max_var |  (Function == "exonic" & IC_w > repel_thresh )  ) ) %>%
  arrange(desc(IC_w)) %>% distinct(gene_id, .keep_all = TRUE)


manh_plot <- ggplot(d_f , aes(x = dist, y = IC_w)) +
  # IC_woints for all variants
  geom_point(data=subset(df_nonexonic,!chr_max_var | IC_w < repel_thresh) , aes(color = as.factor(CHR)), size=0.5*dpi/72, shape=16,alpha=1, fill="blue", stroke=0.5) + #size= 0.2 * dpi/72,
  geom_point(data=subset(df_nonexonic,chr_max_var & IC_w >= repel_thresh) , color="black", shape=23,alpha=1,size= 1*dpi/72, fill="blue", stroke=0.5) +
  #scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 16))  +
  #scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6))  +
  scale_color_manual(values = rep(c("gray20", "gray50"), length.out = length(unique(d_f$CHR)))) +
  new_scale_color() +
#  geom_point(data=df_exonic, aes(shape = factor(chr_max_var),size=factor(chr_max_var), color=factor(chr_max_var)  ),  alpha = 0.6, fill="green", stroke=0.5) +
  geom_point(data=df_exonic, aes(size=factor(chr_max_var) ), color="black", shape=24, alpha = 1, fill="green", stroke=0.5) +
#  scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 17))  +
  scale_size_manual(values = c("TRUE" = 2*dpi/72, "FALSE" = 1 * dpi/72))  +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "green"))  +
  geom_hline(yintercept = y_min, color = "black", linewidth=0.1 * dpi/72) +
  new_scale_color() + new_scale('shape') + new_scale('size') + 
  geom_text_repel(data=d_f_unique_gene_id,
    aes(x = dist, y = IC_w, label = gene_id, color=factor(var_type) ),
      size = 2 * dpi/72, angle=90 , #arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "last"),
#    box.padding = 0,
    max.overlaps = Inf, na.rm = FALSE, parse=TRUE,
    min.segment.length = 0, segment.color = "red",
#   force = 1, #force_pull= , 
    nudge_y = 1,            
#    nudge_x = 300000,           
    direction = "both",
    point.padding = 0, family="Arial") + 
    scale_color_manual(values = c("nonexonic" = "black", "exonic" = "ForestGreen"))+
  scale_x_continuous(breaks = chr_breaks, labels = chr_labels) + coord_cartesian(clip = "off", ylim=c(y_min, y_max)) +
#  scale_shape_manual(values = c("exonic" = 17, "non-exonic" = 16)) +
#  scale_size_manual(values = c("exonic" = 5, "non-exonic" = 1)) +
  labs(x = "Chromosome", y = "IC z-score" ) + #, title = plot_title) +
  theme_minimal() +
  theme(plot.margin = margin(t = 10 , r = 10, b = 10, l = 10, unit="mm"), legend.position = "none",  panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank() ,axis.text.x = element_blank(), plot.background = element_rect(fill = "white"),
text = element_text(family = "Arial", size=3* dpi/72),
    axis.title = element_text(family = "Arial", size=3 * dpi/72),
    axis.text = element_text(family = "Arial", size=3* dpi/72),
    plot.title = element_text(family = "Arial", size=4* dpi/72),
    plot.subtitle = element_text(family = "Arial", size=3* dpi/72),
axis.line.y = element_line(color = "black", linewidth=0.1 * dpi/72),
axis.ticks.y = element_line(color = "black"),
axis.ticks.length.y = unit(0.15, "cm")) +
  geom_text_repel(data = data.frame(x = chr_breaks, y = y_min - 0.05, label = chr_labels),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, na.rm = FALSE,
    direction = "x",
    min.segment.length = 0.2*dpi/72,
    nudge_y = 0,  
#    nudge_x = 0,
    size = 1 * dpi/72,
    angle = 0, family="Arial", segment.color = "gray70",
    box.padding = 0.2, max.time=20,vjust=1,
    point.padding = 0.2,
    max.overlaps = Inf) + 
  scale_y_continuous(labels = label_number(accuracy = 0.1), breaks = seq(0, ceiling(y_max / 0.5) * 0.5 , by = 0.5))
  # limits = c(
#		     floor(min(y_max * sgn, y_min)),
#		     ceiling(max(y_max * sgn, y_min * sgn)) ),
#  geom_hline(yintercept = 7.30103 * sgn , linetype = "dashed", color = "red", linewidth=0.1 * dpi/72) +	
#  geom_hline(yintercept = y_inter*sgn , linetype = "solid", color = "red", linewidth=0.1 * dpi/72) +


#if (contrast!='longitudinal') {plots[[contrast]]  <- plots[[contrast]] + geom_hline(yintercept = y_inter*sgn , linetype = "solid", color = "red", linewidth=0.1 * dpi/72) }
#plots[[contrast]]  <- plots[[contrast]] + geom_hline(yintercept = y_inter*sgn , linetype = "solid", color = "red", linewidth=0.1 * dpi/72)

#ggsave(paste0('manhattan_ic_',ic,'.png'), plot = m_plot, width = 180*4, height = 170*4, units = "mm", dpi = dpi)


df <- data.frame(
  expected = -log10((1:n)/(n+1)),
  observed = -log10(sort(m_p))
)

zzplot <- ggplot(df, aes(expected, observed)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, color="red") +
  labs(x = "Expected -log10(p)",
       y = "Observed -log10(p)") +
  theme_classic()
	
#ggsave(paste0('qqplot_ic_',ic,'.png'), plot = zzplot, width = 180, height = 180, units = "mm", dpi = dpi)


df <- data.frame( z = m_z)
histo <- ggplot(df, aes(x = z)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 100,
                 color = "black",
                 fill = "grey80") +
  stat_function(fun = dnorm,
                args = list(mean = 0, sd = 1),
                color = "blue",
                linetype = "dashed",
                linewidth = 1) +
  labs(x = "z", y = "Density") +
  theme_classic()

#ggsave(paste0('histo_ic_',ic,'.png'), plot = histo, width = 180, height = 180, units = "mm", dpi = dpi)



library(patchwork)
combined <- (manh_plot | (zzplot / histo))  + plot_layout(widths = c(4,1))
ggsave(paste0("manhattan_ic_",ic,".png"), combined, width = 50*25, height = 180*4, units = "mm", dpi = dpi)

print("PNG saved")


}




