#############################################
# spatial_DG_lifespan project
# BayesSpace Spatial Clustering for only k=10
# Anthony Ramnauth, March 28 2023
#############################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(BayesSpace)
    library(ggplot2)
    library(Polychrome)
})

# Load SPE
spe <-
    readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <-
    list(platform = "Visium", is.enhanced = FALSE)

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
colData(spe)$row <-
    colData(spe)$array_row + auto_offset_row[spe$sample_id]
colData(spe)$col <- colData(spe)$array_col

# Run BayesSpace
message("Running spatialCluster()")
Sys.time()
set.seed(12345)
spe <-
    spatialCluster(
        spe,
        use.dimred = "HARMONY",
        q = 10,
        platform = "Visium",
        nrep = 50000
    )
Sys.time()

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", 10)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data", "k10_clustering_results"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()