###################################################################################################
# spatial_DG_lifespan project
# Check spatial coordinates & orientation of individual capture areas within spe object are correct
# Anthony Ramnauth, Apr 11 2022
###################################################################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

library(SpatialExperiment)
library(ggspavis)
library(sessioninfo)
library(here)

# Create directory for QC plots
dir_rdata <- here::here("plots", "QC_plots")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

# Load spe
spe <- readRDS(here::here("processed-data", "02_build_spe","spe.rds"))

# Plot spatial coordinates of within tissue spots
pdf(file = here::here("plots", "QC_plots", "spe_Spatial_Coordinates.pdf"))
plotSpots(spe = spe[, colData(spe)$sample_id == "Br8686"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue")
plotSpots(spe = spe[, colData(spe)$sample_id == "Br2706"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue")
plotSpots(spe = spe[, colData(spe)$sample_id == "Br3942"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue")
plotSpots(spe = spe[, colData(spe)$sample_id == "Br6023"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue")

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
