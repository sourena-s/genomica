# module load 2025 ; module load R/4.4.2-gfbf-2025a

library(corrplot)

dims <- c(5, 10, 15, 20)

#all_list <- lapply(dims, function(dim) {
#  mix <- read.table(
#    paste0("matrix_4d_r2psy_ICA_dim_", dim, "/melodic_mix"),
#    header = FALSE
#  )

#  colnames(mix) <- paste0("R2PSY_IC", seq_len(dim), "_dim", dim)
#  mix
#})

#r2psy_concatenated_matrix <- do.call(cbind, all_list)

#m<-cor(r2psy_concatenated_matrix)

df <- read.table("corr_psy2psy.2", header = FALSE)
colnames(df) <- c("ICNAME", "X", "Y", "val")

df$X <- as.numeric(df$X)
df$Y <- as.numeric(df$Y)

df$X <- factor(df$X, levels = sort(unique(df$X)))
df$Y <- factor(df$Y, levels = sort(unique(df$Y)))

mat <- xtabs(val ~ X + Y, data = df)

df_cell_name <- read.table("corr_psy2psy.2", header = FALSE, nrows=50)[,c(1,3)]
colnames(df_cell_name) <-c("IC","ID")
df_cell_name$IC <-substring(df_cell_name$IC, 12)

colnames(mat) <- df_cell_name$IC
rownames(mat) <- df_cell_name$IC



#all_list <- lapply(dims, function(dim) {
#  mix <- read.table(
#    paste0("matrix_4d_manual_rpoint3_ICA_dim_", dim, "/melodic_mix"),
#    header = FALSE
#  )
#
#  colnames(mix) <- paste0("MANUAL_IC", seq_len(dim), "_dim", dim)
#  mix
#})

#manual_concatenated_matrix <- do.call(cbind, all_list)

#m<-cor(manual_concatenated_matrix)
#corrplot(m,order='hclust')


#m_2_m <- cor(r2psy_concatenated_matrix, manual_concatenated_matrix)

svg("corrplot.svg", width = 12, height = 12, bg = "white")

corrplot(mat,order='hclust',outline=T, addrect=14,rect.col="Red" ,tl.cex = 0.7 )
#corrplot(m_2_m,order='hclust')
dev.off()


