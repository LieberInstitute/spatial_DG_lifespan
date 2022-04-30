###############################
# spatial_DG_lifespan project
# Batch correction
# Anthony Ramnauth, Apr 21 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(harmony)
    library(ggplot2)
    library(ggnewscale)
    library(PCAtools)
    library(spatialLIBD)
    library(sessioninfo)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "batch_correction")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

# Feature selection
dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
head(top_hvgs)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
head(top_hvgs)

# Dimensionality reduction
set.seed(12345)
spe <- runPCA(spe, subset_row = top_hvgs, ncomponents = 50)
set.seed(12345)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")

# Plot UMAP before harmony
pdf(file=here::here("plots", "batch_correction", "DG_UMAP_spe.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id), alpha = 0.5)) +
  geom_point() +
  labs(color = "Capture Area") +
  theme_bw()

dev.off()

# Clustering
set.seed(12345)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)

# Plot graph-based clustering
df <- cbind.data.frame(colData(spe), spatialCoords(spe),
                       reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"))

pdf(file = here::here("plots", "batch_correction", "DG_graph_clusters_spe.pdf"))
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) +
  geom_point(size = 0.2, alpha = 0.5) +
  ggtitle("Graph-based clustering: no batch correction") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_bw() +
  theme(panel.grid = element_blank())

dev.off()

# Plot pre-harmony integrated graph-based clusters onto tissue
pdf(file = here::here("plots", "batch_correction", "DG_tissue_plot_graph_clusters_spe.pdf"))
  ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                     color = label)) +
    facet_wrap(~ sample_id) +
    geom_point(size = 0.1) +
    scale_y_reverse() +
    ggtitle("graph-based clusters pre-harmony") +
    labs(color = "cluster") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

dev.off()


# run Harmony on PCA dimensions to integrate sample IDs

pca_matrix <- reducedDim(spe, "PCA")
sample_ids <- colData(spe)$sample_id
stopifnot(nrow(pca_matrix) == length(sample_ids))

set.seed(123)
harmony_embeddings <- HarmonyMatrix(
  pca_matrix,
  meta_data = sample_ids,
  do_pca = FALSE
)

colnames(harmony_embeddings) <- paste0("HARMONY", seq_len(ncol(harmony_embeddings)))

dim(harmony_embeddings)
head(harmony_embeddings, 2)

# store in SPE object
reducedDims(spe) <- list(
  PCA = reducedDim(spe, "PCA"),
  UMAP = reducedDim(spe, "UMAP"),
  HARMONY = harmony_embeddings
)

# Run UMAP on harmony integrated spe object
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

# Plot UMAP after harmony
pdf(file=here::here("plots", "batch_correction", "DG_UMAP_harmony.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id), alpha = 0.5)) +
  geom_point() +
  labs(color = "Capture Area") +
  theme_bw()

dev.off()

# Clustering on batch corrected spe object
set.seed(12345)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "HARMONY")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)

# Plot graph-based clustering
df <- cbind.data.frame(colData(spe), spatialCoords(spe),
                       reducedDim(spe, "HARMONY"), reducedDim(spe, "UMAP.HARMONY"))

pdf(file=here::here("plots", "batch_correction", "DG_graph_clusters_harmony.pdf"))
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) +
  geom_point(size = 0.2, alpha = 0.5) +
  ggtitle("Graph-based clustering: batch corrected") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_bw() +
  theme(panel.grid = element_blank())

dev.off()

# Plot post-harmony integrated graph-based clusters onto tissue
pdf(file = here::here("plots", "batch_correction", "DG_tissue_plot_graph_clusters_harmony.pdf"))
  ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                     color = label)) +
    facet_wrap(~ sample_id) +
    geom_point(size = 0.1) +
    scale_y_reverse() +
    ggtitle("graph-based clusters post-harmony") +
    labs(color = "cluster") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

dev.off()

# Create directory for harmony integrated spe
dir_rdata <- here::here("processed-data", "harmony_processed_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

saveRDS(spe, file = here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
