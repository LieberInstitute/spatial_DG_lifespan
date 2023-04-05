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

# Set gene names as row names for common gene ID with snRNAseq data
rownames(spe) <- rowData(spe)$gene_name

# Create martix and df for createCARDObject arguments

spatial_cou <- assays(spe)$counts

stopifnot(colnames(spatial_cou) == rownames(colData(spe)))
stopifnot(spe$key == colData(spe)$key)
colnames(spatial_cou) <- spe$key

# The spatial locations for each sample might need to be offset for proper spatial auto-correlation
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
colData(spe)$row <-
    colData(spe)$array_row + auto_offset_row[spe$sample_id]
colData(spe)$col <- colData(spe)$array_col

spatial_loc <- data.frame(
    x = colData(spe)$col,
    y = colData(spe)$row,
    row.names = spe$key
)

stopifnot(colnames(spatial_cou) == rownames(spatial_loc))

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

saveRDS(CARD_obj, file = here::here("processed-data", "Cell_Type_Deconvolution", "CARD_obj_sestan.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
