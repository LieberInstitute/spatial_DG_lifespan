###############################
# spatial_DG_lifespan project
# Plots for Manual Annotations
# Anthony Ramnauth, Oct 07 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(RColorBrewer)
    library(ggplot2)
    library(Polychrome)
    library(ggspavis)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# order spe observations according to age
spe <- spe[, order(spe$age)]

# Assigning names to the colors
man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL_CA1" = "#16FF32", "PCL_CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD",
    "SUB" = "#AA0DFE", "THAL" = "navy")

# Plot Manual Annotations onto tissue
vis_grid_clus(
    spe = spe,
    clustervar = "ManualAnnotation",
    pdf = here(
        "plots",
        "manual_annotations",
        "ManualAnnotations.pdf"),
    sort_clust = FALSE,
    colors = man_colors,
    spatial = TRUE,
    point_size = 2,
    image_id = "lowres",
    alpha = 0.5)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
