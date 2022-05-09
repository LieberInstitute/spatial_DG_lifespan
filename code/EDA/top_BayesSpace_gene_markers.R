##########################################
# spatial_DG_lifespan project
# Top Marker Genes for BayesSpace Clusters
# Anthony Ramnauth, May 09 2022
##########################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(SpatialExperiment)
    library(BayesSpace)
    library(scran)
    library(scater)
    library(pheatmap)
})

# Create directory for Top Marker Genes plots
dir_plots <- here::here("plots", "top_BayesSpace_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "clustering_results"),
  prefix = ""
)

spe$bayesSpace_harmony_8 <- as.factor(spe$bayesSpace_harmony_8)

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

# Test for marker genes
markers <- findMarkers(spe, groups = spe$bayesSpace_harmony_8, test = "binom", direction = "up")

# Returns a list with one DataFrame per cluster
markers

# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 3 since BayesSpace tissue plot looks like GCL
interesting_3 <- markers[[3]]
best_set_3 <- interesting_3[interesting_3$Top <= 5, ]
logFCs_3 <- getMarkerEffects(best_set_3)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster3_GCL.pdf"))
pheatmap(logFCs_3, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 3 (GCL)")
dev.off()

# Plot log-transformed normalized expression of top genes for GCL cluster
top_genes_3 <- head(rownames(interesting_3))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster3_GCL.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_3, colour_by = "bayesSpace_harmony_8")
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 4 since BayesSpace tissue plot looks like SGZ
interesting_4 <- markers[[4]]
best_set_4 <- interesting_4[interesting_4$Top <= 5, ]
logFCs_4 <- getMarkerEffects(best_set_4)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster4_SGZ.pdf"))
pheatmap(logFCs_4, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 4 (SGZ)")
dev.off()

# Plot log-transformed normalized expression of top genes for SGZ cluster
top_genes_4 <- head(rownames(interesting_4))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster4_SGZ.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_4, colour_by = "bayesSpace_harmony_8")
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 5 since BayesSpace tissue plot looks like CA4
interesting_5 <- markers[[5]]
best_set_5 <- interesting_5[interesting_5$Top <= 5, ]
logFCs_5 <- getMarkerEffects(best_set_5)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster5_CA4.pdf"))
pheatmap(logFCs_5, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 5 (CA4)")
dev.off()

# Plot log-transformed normalized expression of top genes for CA4 cluster
top_genes_5 <- head(rownames(interesting_5))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster5_CA4.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_5, colour_by = "bayesSpace_harmony_8")
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 7 since BayesSpace tissue plot looks like ML
interesting_7 <- markers[[7]]
best_set_7 <- interesting_7[interesting_7$Top <= 5, ]
logFCs_7 <- getMarkerEffects(best_set_7)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster7_ML.pdf"))
pheatmap(logFCs_7, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 7 (ML)")
dev.off()

# Plot log-transformed normalized expression of top genes for ML cluster
top_genes_7 <- head(rownames(interesting_7))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster7_ML.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_7, colour_by = "bayesSpace_harmony_8")
dev.off()

# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 6 since BayesSpace tissue plot looks like CA3
interesting_6 <- markers[[6]]
best_set_6 <- interesting_6[interesting_6$Top <= 5, ]
logFCs_6 <- getMarkerEffects(best_set_6)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster6_CA3.pdf"))
pheatmap(logFCs_6, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 6 (CA2-3)")
dev.off()

# Plot log-transformed normalized expression of top genes for CA3 cluster
top_genes_6 <- head(rownames(interesting_6))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster6_CA3.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_6, colour_by = "bayesSpace_harmony_8")
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
