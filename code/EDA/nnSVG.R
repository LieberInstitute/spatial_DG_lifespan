################################################################
# spatial_DG_lifespan project
# nnSVG per Capture Area, Average ranks, & BayesSpace covariates
# Anthony Ramnauth, July 11 2022
################################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(scater)
    library(scran)
    library(nnSVG)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(sessioninfo)
})

spe <- readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

# Create vector of samples for nnSVG on whole tissue
sample_ids <- c(
    "Br1412",
    "Br2706",
    "Br2720",
    "Br3942",
    "Br5242",
    "Br5699_new",
    "Br6023",
    "Br6129_new",
    "Br6299_new",
    "Br6522",
    "Br8181",
    "Br8195",
    "Br8533",
    "Br8667",
    "Br8686",
    "Br8700"
)

# Run nnSVG once per sample whole tissue and store lists of top SVGs

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {

    # select sample_id
    ix <- colData(spe)$sample_id == sample_ids[s]
    spe_sub <- spe[, ix]

    # run nnSVG filtering for mitochondrial gene and low-expressed genes
    spe_sub <- filter_genes(spe_sub)

	# re-calculate library size factors
	spe_sub <- computeLibraryFactors(spe_sub)

    # re-calculate logcounts after filtering
    spe_sub <- logNormCounts(spe_sub)

    # run nnSVG
    set.seed(12345)
    spe_sub <- nnSVG(spe_sub, n_threads = 10)

    # store whole tissue results
    res_list[[s]] <- rowData(spe_sub)
}

# directory to save whole tissue results
dir_outputs <- here("processed-data", "nnSVG", "whole_tissue")

# save whole tissue nnSVG results
fn_out <- file.path(dir_outputs, "DG_nnSVG_results")
saveRDS(res_list, paste0(fn_out, ".rds"))
save(res_list, file = paste0(fn_out, ".RData"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
