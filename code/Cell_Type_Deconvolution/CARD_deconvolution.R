###################################
# spatial_DG_lifespan project
# Cell-Type Deconvolution with CARD
# Anthony Ramnauth, Oct 14 2022
###################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

rm(list=ls())

suppressPackageStartupMessages({
    library(SingleCellExperiment)
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
sce <- load(file = here::here("processed-data", "SCE_HPC-n3_tran_etal.rda"))

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

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
    ct.varname = "cellType",
    ct.select = unique(sce$cellType),
    sample.varname = "donor",
    minCountGene = 100,
    minCountSpot = 5)

# Run CARD deconvolution

CARD_obj <- CARD_deconvolution(CARD_obj)

saveRDS(CARD_deconvo, file = here::here("CARD_obj.rds"))

dim(CARD_deconvo@Proportion_CARD)
dim(spe)

setdiff(rownames(spatial_loc), rownames(CARD_deconvo@Proportion_CARD))

# There are 9 Visium spots missing, fill in those blanks!

emptyNaDF <- data.frame(matrix(NA,nrow = 9, ncol = 23))
rownames(emptyNaDF) <- setdiff(rownames(spatial_loc), rownames(CARD_deconvo@Proportion_CARD))
colnames(emptyNaDF) <- colnames(CARD_deconvo@Proportion_CARD)

cell_props <- rbind(CARD_deconvo@Proportion_CARD, emptyNaDF)

stopifnot(rownames(cell_props) == rownames(spatial_loc))

# Reorder dataframe to match the order of spatial barcodes

cell_props <-  cbind(cell_props, rownames(cell_props))
cell_props <- rename(cell_props, "Barcodes" = "rownames(blah)")
spatial_loc2 <- cbind(spatial_loc, rownames(spatial_loc))
spatial_loc2 <- rename(spatial_loc2, "Barcodes" = "rownames(spatial_loc)")

cell_props <- cell_props[order(match(cell_props$Barcodes, spatial_loc2$Barcodes)), ]
cell_props$Barcodes <- NULL

stopifnot(rownames(cell_props) == rownames(spatial_loc))
cell_props[is.na(cell_props)] = 0

colData(spe) <- cbind(colData(spe), cell_props)

saveRDS(spe, file = here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
