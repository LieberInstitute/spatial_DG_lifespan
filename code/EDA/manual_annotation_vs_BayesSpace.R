###########################################
# spatial_DG_lifespan project
# Comparing Manual Annotation vs BayesSpace
# Anthony Ramnauth, Oct 21 2022
###########################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(ggplot2)
    library(ggnewscale)
    library(spatialLIBD)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL_CA1" = "#16FF32", "PCL_CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD", "SUB" = "#AA0DFE")

bay_colors <- c("1" = "#E4E1E3", "2" = "#FEAF16", "3" = "#FE00FA", "4" = "#B00068",
    "5" = "#DEA0FD", "6" = "#5A5156", "7" = "#3283FE", "8" = "#1CFFCE")

# Plot UMAP after harmony
pdf(file = here::here("plots", "manual_annotations", "MA_vs_BS_UMAP_harmony.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$bayesSpace_harmony_8), alpha = 0.01)) +
    geom_point() +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace Cluster") +
    theme_bw()

ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$ManualAnnotation), alpha = 0.01)) +
    geom_point() +
    scale_color_manual(values = man_colors) +
    labs(color = "Manual Annotation") +
    theme_bw()

dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
