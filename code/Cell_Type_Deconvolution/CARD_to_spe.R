###################################
# spatial_DG_lifespan project
# CARD to spe
# Anthony Ramnauth, Dec 19 2022
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

# Remove CP cluster
spe = spe[, which(spe$bayesSpace_harmony_8 != "5")]

CARD_obj <- readRDS(file = here::here("processed-data", "Cell_Type_Deconvolution", "CARD_obj_sestan.rds"))

dim(CARD_obj@Proportion_CARD)
dim(spe)

# Create martix and df for createCARDObject arguments

spatial_cou <- assays(spe)$counts

spatial_loc <- spatialCoords(spe)

stopifnot(colnames(spatial_cou) == rownames(spatial_loc))

spatial_loc <- data.frame(spatial_loc)
colnames(spatial_cou) <- rownames(spatial_loc)

colnames(spatial_loc) <- c("x", "y")

setdiff(rownames(spatial_loc), rownames(CARD_obj@Proportion_CARD))

# There are 5 Visium spots missing, fill in those blanks!

emptyNaDF <- data.frame(matrix(NA,nrow = 5, ncol = 31))
rownames(emptyNaDF) <- setdiff(rownames(spatial_loc), rownames(CARD_obj@Proportion_CARD))
colnames(emptyNaDF) <- colnames(CARD_obj@Proportion_CARD)

cell_props <- rbind(CARD_obj@Proportion_CARD, emptyNaDF)

cell_props <-  cbind(cell_props, rownames(cell_props))
cell_props <- rename(cell_props, "Barcodes" = "rownames(cell_props)")
spatial_loc2 <- cbind(spatial_loc, rownames(spatial_loc))
spatial_loc2 <- rename(spatial_loc2, "Barcodes" = "rownames(spatial_loc)")

cell_props <- cell_props[order(match(cell_props$Barcodes, spatial_loc2$Barcodes)), ]
cell_props$Barcodes <- NULL

stopifnot(rownames(cell_props) == rownames(spatial_loc))
cell_props[is.na(cell_props)] = 0

colData(spe) <- cbind(colData(spe), cell_props)

# Find dominant celltype for each spot
p <- cell_props %>%
 rowwise() %>%
 mutate(row_max = names(.)[which.max(c_across(everything()))])

colData(spe)$dominant_cell_types <- p$row_max

saveRDS(spe, file = here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_CARD.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
