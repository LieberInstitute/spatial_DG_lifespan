###########################################
# spatial_DG_lifespan project
# Comparing Manual Annotation vs BayesSpace
# Anthony Ramnauth, Oct 21 2022
###########################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(ggplot2)
    library(ggnewscale)
    library(spatialLIBD)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL_CA1" = "#16FF32", "PCL_CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD",
    "SUB" = "#AA0DFE", "THAL" = "navy")

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

# Plot UMAP after harmony
pdf(file = here::here("plots", "manual_annotations", "MA_vs_BS_UMAP_harmony.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$bayesSpace_harmony_10), alpha = 0.01)) +
    geom_point() +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace Cluster") +
    theme_bw() +
    theme_classic()

ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$ManualAnnotation), alpha = 0.01)) +
    geom_point() +
    scale_color_manual(values = man_colors) +
    labs(color = "Manual Annotation") +
    theme_bw() +
    theme_classic()

dev.off()
