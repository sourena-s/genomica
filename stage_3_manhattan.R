library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggnewscale)
library(svglite)
library(patchwork)

#stat_file= commandArgs(trailingOnly = T)[1]
#output_png= commandArgs(trailingOnly = T)[2]

plots<-list()
plot_title <- "Genomic IC"

repel_pthresh <- 8# ggrepel P-value threshold

y_max <- 15
y_min <- 0 
y_inter <- -log10((5)*10^-8)
exonic_pval_threshold <- -log10((5)*10^-8)
top_chr_var_plot_pthresh <- 8
  

stat_file=paste0('/data/clusterfs/lag/users/sousoh/impact_modes/ica/', contrast, '/merged/main-alldim-minp.fil')

annovar_cols <- c("Function", "Gene", "CHR", "Start", "End", "Ref", "Alt", "SNP", "VariantTag", "P_annovar")
annovar1 <- fread( paste0('sig_', contrast_str, '.annovar.variant_function.txt') , header = FALSE, sep = " ")
colnames(annovar1) <- annovar_cols
annovar1 <- annovar1[annovar1$Function != "exonic",]
annovar2 <-  fread( paste0('sig_', contrast_str, '.annovar.exonic_variant_function.txt') , header = FALSE, sep = " ")
colnames(annovar2) <- annovar_cols
annovar <-rbind(annovar1, annovar2)

#intronic        CPED1   7 120927437 120927437 A AC rs112901683 7:120927437[imp]A,AC 57.4682
#ID CHROM POS ORIG_RSID REF ALT A1 NEG_LOG10_P
#1:756604[call]A,G 1 756604 rs3131962 A G A 2.98305

print(stat_file)


a<-as.data.frame(fread(stat_file,header=T))
b<-a[!is.na(a$NEG_LOG10_P) & a$NEG_LOG10_P >= y_min,]


pthresh_sparse <- -log10((5)*10^-8)

nbins <- 20
n_per_bin <- 200
print(dim(b))
b <- bind_rows(b %>% filter(NEG_LOG10_P > pthresh_sparse), b %>% filter(NEG_LOG10_P > 0 & NEG_LOG10_P <= pthresh_sparse) %>% mutate(bin = cut(NEG_LOG10_P, breaks = nbins, include.lowest = TRUE)) %>% group_by(CHROM,bin) %>% filter(n() > 0) %>% slice_sample(n = n_per_bin, replace = TRUE) %>% ungroup()) %>% select(-bin) ## %>% arrange(NEG_LOG10_P)

print(dim(b))
colnames(b) <- c("VariantTag","CHR","BP","SNP","A1","A2","A3","P")
######DEBUG
qc_1<-as.data.frame(fread(qc_file,header=F,sep=" "))[,c(2,5)]
qc_2<-as.data.frame(fread(qc_vartag,header=F,sep=" "))[,1]
qc<-cbind(qc_1,qc_2)
colnames(qc) <- c("SNP_QC","QC","VariantTag")
print('Joining GWAS and QC list')
m_joined <- left_join(b,qc,by="VariantTag")


print('Merging GWAS and annotations')

m_joined<- merge(m_joined, annovar[, c("VariantTag", "Gene", "Function")], by="VariantTag", all.x = TRUE)

m_joined$Gene_OLD <- m_joined$Gene
m_joined$Gene <- sub("\\(.*", "", m_joined$Gene_OLD)

top_exonic_snps_per_gene <- m_joined[m_joined$P > exonic_pval_threshold  ,] %>%
  filter(!is.na(Gene) & (Function == "nonsynonymous" | Function == "frameshift" | Function == "nonframeshift" | Function == "stopgain") ) %>%
  group_by(Gene) %>%
  slice_max(order_by = P, n = 1, with_ties = FALSE) %>% select(SNP)

print(paste('Found exonic variants', top_exonic_snps_per_gene))
 
top_snps_per_gene <- m_joined[m_joined$P > repel_pthresh,] %>%
  filter(!is.na(Gene) & Function != "intergenic") %>%
  group_by(Gene) %>%
  slice_max(order_by = P, n = 1, with_ties = FALSE) %>% select(SNP) %>% filter(!SNP %in% top_exonic_snps_per_gene)


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
gap <- 40000000
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
  chr_subset <- which(d_f$CHR == i & !is.na(d_f$P) & d_f$P > top_chr_var_plot_pthresh )
  if (length(chr_subset) > 0) {
    top_idx <- chr_subset[which.max(d_f$P[chr_subset])]
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

if (grepl('baseline', contrast_str) || contrast_str == 'diff') {
  rows_to_update <- which(d_f$gene_id == 'rs1080066' & !is.na(d_f$gene_id))
  if (length(rows_to_update) > 0) {
    d_f$gene_id[rows_to_update] <- "rs1080066 ~'(' * italic('THBS1') * ')'"
  }
}

max_length <- max(nchar(d_f$gene_id), na.rm = TRUE)


label_y_position <- -0.1 + y_min

d_f$CHR <- as.numeric(d_f$CHR)

dpi=300 ##use 72 for SVG

d_f$gene_id <- gsub("_", ":", d_f$gene_id)

d_f <- subset(d_f, P>0 & is.finite(dist) & is.finite(P) ) 

d_f$var_type <- ifelse(d_f$Function %in% c("nonsynonymous", "frameshift", "nonframeshift", "stopgain"), "exonic", "nonexonic")

df_nonexonic <- subset(d_f,  ( Function != "nonsynonymous" & Function != "frameshift" & Function != "nonframeshift" & Function != "stopgain" ) | P <= exonic_pval_threshold) 
df_exonic    <- subset(d_f,  P > exonic_pval_threshold & !is.na(Function) & (Function == "nonsynonymous" | Function == "frameshift" | Function == "nonframeshift" | Function == "stopgain") |
		                 (chr_max_var &  (Function == "nonsynonymous" | Function == "frameshift" | Function == "nonframeshift" | Function == "stopgain")  )
		       ) 

chr_breaks <- chr_breaks[c(1,3,5,7,9,11,13,15,17,19,21,23,24,25)]
chr_labels <- chr_labels[c(1,3,5,7,9,11,13,15,17,19,21,23,24,25)]

d_f_unique_gene_id <- d_f %>%
  filter( (P > repel_pthresh & !is.na(gene_id)) |
 ( chr_max_var &  (Function == "nonsynonymous" | Function == "frameshift" | Function == "nonframeshift" | Function == "stopgain") ) ) %>%
  arrange(desc(P)) %>% distinct(gene_id, .keep_all = TRUE)

sgn = 1 #inverted Manhattan plot: -1

m_plot <- ggplot(d_f , aes(x = dist, y = P*sgn)) +
  # Points for all variants
  geom_point(data=subset(df_nonexonic,!chr_max_var | P < repel_pthresh) , aes(color = as.factor(CHR),size=factor(P > -log10((5/848)*10^-8) ) ), shape=16,alpha=0.6, fill="blue", stroke=0.5) + #size= 0.2 * dpi/72,
  geom_point(data=subset(df_nonexonic,chr_max_var & P >= repel_pthresh) , color="black", shape=23,alpha=1,size= 0.6 * dpi/72, fill="blue", stroke=0.5) +
  #scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 16))  +
  #scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6))  +
  scale_size_manual(values = c("TRUE" = 0.2* dpi/72, "FALSE" = 0.2 * dpi/72))  +
  scale_color_manual(values = rep(c("gray20", "gray50"), length.out = length(unique(d_f$CHR)))) +
  new_scale('size') + 
  new_scale_color() +
#  geom_point(data=df_exonic, aes(shape = factor(chr_max_var),size=factor(chr_max_var), color=factor(chr_max_var)  ),  alpha = 0.6, fill="green", stroke=0.5) +
  geom_point(data=df_exonic, aes(size=factor(chr_max_var) ), color="black", shape=24.5 - sgn*0.5, alpha = 1, fill="green", stroke=0.5) +
#  scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 17))  +
  scale_size_manual(values = c("TRUE" = 0.5*dpi/72, "FALSE" = 0.4 * dpi/72))  +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "green"))  +
  geom_hline(yintercept = y_min, color = "black", linewidth=0.1 * dpi/72) +
  new_scale_color() + new_scale('shape') + new_scale('size') + 
  geom_text_repel(data=d_f_unique_gene_id,
    aes(x = dist, y = P*sgn, label = gene_id, color=factor(var_type) ),
      size = 1 * dpi/72, angle=ifelse(contrast=='longitudinal', 90, 90) , #arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "last"),
#    box.padding = 0,
    max.overlaps = Inf, na.rm = FALSE, parse=TRUE,
    min.segment.length = 0, segment.color = "red",
#   force = 1, #force_pull= , 
    nudge_y = ifelse(contrast=='longitudinal', 1*sgn, 0.5 * sgn),            
#    nudge_x = 300000,           
    direction = "both",
    point.padding = 0, family="Arial"
  ) + scale_color_manual(values = c("nonexonic" = "black", "exonic" = "ForestGreen"))+
  scale_x_continuous(breaks = chr_breaks, labels = chr_labels) + coord_cartesian(clip = "off", ylim=c(min( y_max*sgn, y_min), max(y_max*sgn, y_min*sgn))) +
#  scale_shape_manual(values = c("exonic" = 17, "non-exonic" = 16)) +
#  scale_size_manual(values = c("exonic" = 5, "non-exonic" = 1)) +
  labs(x = "Chromosome", y = expression(-Log(italic(P))) ) + #, title = plot_title) +
  theme_minimal() +
  theme(plot.margin = margin(t = 10 , r = 10, b = 10, l = 10, unit="mm"), legend.position = "none",  panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank() ,axis.text.x = element_blank(), plot.background = element_rect(fill = "white"),
text = element_text(family = "Arial", size=3* dpi/72),
    axis.title = element_text(family = "Arial", size=3 * dpi/72),
    axis.text = element_text(family = "Arial", size=3* dpi/72),
    plot.title = element_text(family = "Arial", size=4* dpi/72),
    plot.subtitle = element_text(family = "Arial", size=3* dpi/72),
axis.line.y = element_line(color = "black", linewidth=0.1 * dpi/72),
axis.ticks.y = element_line(color = "black"),
axis.ticks.length.y = unit(0.15, "cm")
	) +
  geom_hline(yintercept = 7.30103 * sgn , linetype = "dashed", color = "red", linewidth=0.1 * dpi/72) +	
  geom_hline(yintercept = y_inter*sgn , linetype = "solid", color = "red", linewidth=0.1 * dpi/72) +
  geom_text_repel(data = data.frame(x = chr_breaks, y = label_y_position, label = chr_labels),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, na.rm = FALSE,
    direction = "x",
    min.segment.length = 0,
    nudge_y = ifelse(contrast=='longitudinal', -15 * (y_max - y_min) /250 , -15*sgn),  
#    nudge_x = 0,
    size = 1 * dpi/72,
    angle = 0, family="Arial", segment.color = "gray70",
    box.padding = 0.2, max.time=20,vjust=1,
    point.padding = 0.2,
    max.overlaps = Inf) + 
  scale_y_continuous(labels = function(x) abs(round(x)),
		     breaks = seq(floor(min(y_max * sgn, 0)),
      	             ceiling(max(y_max * sgn, 0 * sgn)), by = ifelse(contrast=='longitudinal', 5 , 50)  ))
  # limits = c(
#		     floor(min(y_max * sgn, y_min)),
#		     ceiling(max(y_max * sgn, y_min * sgn)) ),


#if (contrast!='longitudinal') {plots[[contrast]]  <- plots[[contrast]] + geom_hline(yintercept = y_inter*sgn , linetype = "solid", color = "red", linewidth=0.1 * dpi/72) }
#plots[[contrast]]  <- plots[[contrast]] + geom_hline(yintercept = y_inter*sgn , linetype = "solid", color = "red", linewidth=0.1 * dpi/72)

ggsave(paste0(contrast_str, '.png'), plot = m_plot, width = 180, height = 170, units = "mm", dpi = dpi)
print("PNG saved")
ggsave(paste0(contrast_str, '.svg'), plot = m_plot, width = 180, height = 170, units = "mm")
print("SVG saved")




