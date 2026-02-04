# install.packages(c("data.table", "uwot", "ggplot2", "ggrepel"))

library(data.table)
library(uwot)
library(ggplot2)
library(ggrepel)

# ---- Files ----
#mix_file   <- "matrix_4d_r2psy_ICA_dim_50/melodic_mix"
#label_file <- "all_h5_rho_point5.pheno_text"

mix_file   <- "matrix_4d_manual_rpoint3_ICA_dim_20/melodic_mix"
label_file <- "curated_phenotype_list_491.pheno_text"

# ---- Read features ----
X <- fread(mix_file, header = FALSE, data.table = FALSE)
X <- as.matrix(X)
storage.mode(X) <- "double"
#stopifnot(ncol(X) == 50)

# ---- Read labels ----
labels <- trimws(readLines(label_file, warn = FALSE))
labels[labels == ""] <- NA
stopifnot(length(labels) == nrow(X))

# ---- UMAP ----
set.seed(1)
X_scaled <- scale(X)

um <- umap(
  X_scaled,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  n_components = 2
)

df <- data.frame(
  UMAP1 = um[, 1],
  UMAP2 = um[, 2],
  label = labels
)

# ---- Select N non-NA labels ----
N <- 100

valid <- which(!is.na(df$label))

score <- abs(scale(df$UMAP1)) + abs(scale(df$UMAP2))
idx <- valid[order(score[valid], decreasing = TRUE)][seq_len(min(N, length(valid)))]

df$label_n <- NA
df$label_n[idx] <- df$label[idx]

# ---- Plot ----
p <- ggplot(df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1) +
  geom_text_repel(
    aes(label = label_n),
    na.rm = TRUE,
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0
  ) +
  theme_classic() +
  labs(title = "UMAP of melodic_mix (50 features)",
       x = "UMAP-1", y = "UMAP-2")

print(p)

