###############################
# spatial_DG_lifespan project
# clustering of sce object
# Anthony Ramnauth, Nov 05 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(bluster)
    library(ggplot2)
    library(viridisLite)
    library(dplyr)
    library(tidyr)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_harmony.rds"))

# Custom clustering algorithm and parameters from both OSCA & Franjic et al. 2021 NEURON
# two-stage clustering algorithm using high-resolution k-means and graph-based clustering
# Will be iterative for each broad cell type clusters

set.seed(12345)
clus <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 2000),
    second = NNGraphParam(k = 25)
  )
)

table(clus)
#clus
#    1     2     3     4     5     6     7     8
#60596 57896 12802 27202 18846  5697 12702 13534

colLabels(sce) <- clus

table(colLabels(sce), colData(sce)$Dataset)
#    Franjic_etal_2022 Zhong_etal_2020 Zhou_etal_2022
#  1             19408            8677          32511
#  2             26494            2515          28887
#  3              4908             222           7672
#  4             12954            1990          12258
#  5              2603            9117           7126
#  6              1863             391           3443
#  7             12267             381             54
#  8              6867            1752           4915

# Plot UMAP after clustering
pdf(file = here::here("plots", "sce_plots", "DG_UMAP_initial_clusters_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label",
    point_alpha = 0.3, point_size = 0.5) +
  ggtitle("Initial broad unsupervised clustering")

dev.off()


# Save new sce object
saveRDS(sce, file = here::here("processed-data", "sce", "sce_clustered.rds"))
