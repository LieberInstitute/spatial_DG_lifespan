###################################
# spatial_DG_lifespan project
# Cell-Type Deconvolution with CARD
# Anthony Ramnauth, Oct 14 2022
###################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

rm(list=ls())

suppressPackageStartupMessages({
    library(SingleCellExperiment)
	library(spatialLIBD)
    library(pbmcapply)
    library(CARD)
    library(here)
    library(SpatialExperiment)
    library(MuSiC)
    library(dplyr)
    library(scatterpie)
    library(sessioninfo)
})


# Load SCE
sce <- readRDS(file = here::here("processed-data", "sce", "sce_sestan_DG_final.rds"))

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

# Remove CP cluster
spe = spe[, which(spe$bayesSpace_harmony_8 != "5")]

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

# Create martix and df for createCARDObject arguments

spatial_cou <- assays(spe)$counts

spatial_loc <- spatialCoords(spe)

stopifnot(colnames(spatial_cou) == rownames(spatial_loc))

spatial_loc <- data.frame(spatial_loc)
colnames(spatial_cou) <- rownames(spatial_loc)

colnames(spatial_loc) <- c("x", "y")

# Create CARD object

CARD_obj <- createCARDObject(
    sc_count = assays(sce)$counts,
    sc_meta = colData(sce),
    spatial_count = spatial_cou,
    spatial_location = spatial_loc,
    ct.varname = "Cell_Type",
    ct.select = unique(sce$Cell_Type),
    sample.varname = "sample_ID",
    minCountGene = 100,
    minCountSpot = 5)

# Run CARD deconvolution

CARD_obj <- CARD_deconvolution(CARD_obj)

saveRDS(CARD_obj, file = here::here("CARD_obj_sestan.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
