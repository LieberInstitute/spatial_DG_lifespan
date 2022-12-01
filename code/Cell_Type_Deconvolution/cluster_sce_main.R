###############################
# spatial_DG_lifespan project
# clustering of sce object
# Anthony Ramnauth, Nov 05 2022
###############################


suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(bluster)
    library(ggplot2)
    library(viridisLite)
    library(ggVennDiagram)
    library(dplyr)
    library(tidyr)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_harmony.rds"))

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and graph-based clustering

set.seed(12345)
clus <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 1000),
    second = NNGraphParam(k = 10, cluster.fun = "leiden")
  )
)

table(clus)

colLabels(sce) <- clus

table(colLabels(sce), colData(sce)$Dataset)

# Plot UMAP after clustering
pdf(file = here::here("sce_plots", "DG_UMAP_clusters_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label",
    point_alpha = 0.3, point_size = 0.5) +
  ggtitle("Unsupervised clustering")

dev.off()


# Save new sce object
saveRDS(sce, file = here::here("sce_objects", "sce_clustered.rds"))
