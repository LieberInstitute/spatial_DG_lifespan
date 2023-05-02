#######################################
# spatial_hpc project
# Saving PRECAST clusters to spe object
# Anthony Ramnauth, Dec 13 2022
#######################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library(here)
})

# Load the PRECAST objects (looks like you can only do one at a time)
load(file = here::here("processed-data", "Seurat", "PRECASTObj_10_SVG.Rdata"))

# Taking this bit from Maddy's script, not sure I need it.
resList <- PRECASTObj@resList
PRECASTObj <- selectModel(PRECASTObj)

# Convert to seurat object for easier format to explore
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
seuInt

# Quickly plot from Seurat object to verify accuracy of spatial clustering
# before adding to spe object
cols <- Polychrome::palette36.colors(15)
names(cols) <- sort(unique(seuInt$cluster))

pdf(file = here::here("plots", "PRECAST_plots", "K15.pdf"), width = 12, height = 28)

SpaPlot(seuInt, batch = NULL, item = "cluster", point_size = 1, cols = cols)

dev.off()

# If clusters look accurate move to adding PRECAST clusters to spe object

# Load SPE
spe <-
    readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

dim(spe)

dim(seuInt@meta.data)

k_15 <- data.frame(
    cluster = seuInt@meta.data$cluster,
    barcodes = rownames(seuInt@meta.data)
)

stopifnot(rownames(k_15) == colnames(spe))

rownames(k_15) <- rownames(seuInt@meta.data)

# Reorder to match colnames(spe)
k_15 <- k_15[order(match(rownames(k_15), colnames(spe))), ]

# Add column of PRECAST cluster assignments to colData(spe)
colData(spe)$PRECAST_k15 <- k_15$cluster

cols <- Polychrome::palette36.colors(10)

# Plot using spatialLIBD to have tissue in background
vis_grid_clus(
    spe = spe,
    clustervar = paste0("bayesSpace_harmony_", k),
    pdf = here(
        "plots",
        "PRECAST_plots",
        "PRECAST_k8.pdf"),
    sort_clust = FALSE,
    colors = cols,
    spatial = TRUE,
    point_size = 2,
    image_id = "lowres",
    alpha = 0.5
)

saveRDS(spe,
    file = here::here("processed-data", "QC_processed_spe", "spe_PRECASTk8.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
