##########################################
# spatial_DG_lifespan project
# Top Marker Genes for BayesSpace Clusters
# Anthony Ramnauth, May 05 2022
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

# set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

# test for marker genes
markers <- findMarkers(spe, groups = spe$bayesSpace_harmony_8, test = "binom", direction = "up")

# returns a list with one DataFrame per cluster
markers

# plot log-fold changes for one cluster over all other clusters
# selecting cluster 6 since BayesSpace tissue plot looks like CA3
interesting <- markers[[6]]
best_set <- interesting[interesting$Top <= 5, ]
logFCs <- getMarkerEffects(best_set)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFCvsCluster6.pdf"))
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101), main = "logFC vs. Cluster 6")
dev.off()

# plot log-transformed normalized expression of top genes for one cluster
top_genes <- head(rownames(interesting), 10)

spe$bayesSpace_harmony_8 <- as.factor(spe$bayesSpace_harmony_8)

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes, colour_by = "bayesSpace_harmony_8")
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
