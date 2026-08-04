library(uwot)
library(ggplot2)
library(plotly)
library(data.table)
library(htmlwidgets)

melodic_mix <- fread('matrix_4d_manual_mergedmask_mergedmask_ICA_dim_20_nl_pow3/melodic_mix',header=F)

labels <- readLines("curated_phenotype_list_491.pheno_text")

um <- umap(melodic_mix, n_components = 2)

df <- data.frame(
  UMAP1 = um[,1],
  UMAP2 = um[,2],
  label = labels
)

p <- ggplot(df, aes(UMAP1, UMAP2, text = label)) +
  geom_point(size = 0.5, alpha = 1) +
  theme_minimal()

interactive_plot <- ggplotly(p, tooltip = "text")  %>% layout(hoverlabel = list(font = list(size = 8)))


htmlwidgets::saveWidget(interactive_plot, "umap_plot.html", selfcontained = FALSE)


ggplot(df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1) +
  geom_text_repel(
    aes(label = label),
    size = 2,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  theme_classic()


