##############################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked across age_bin
# Anthony Ramnauth, Sept 21 2022
##############################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(sessioninfo)
    library(SingleCellExperiment)
    library(rafalib)
    library(limma)
    library(edgeR)
    library(scran)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Format spe object for DE models

colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
    "age",
    "age_bin",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells"
))]

colData(spe_pseudo)

colData(spe_pseudo)$ncells <- as.numeric(colData(spe_pseudo)$ncells)
colData(spe_pseudo)$race <- as.factor(colData(spe_pseudo)$race)
colData(spe_pseudo)$sample_id <- as.factor(colData(spe_pseudo)$sample_id)
colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)

colData(spe_pseudo)

# Drop things we don't need
spatialCoords(spe_pseudo) <- NULL
imgData(spe_pseudo) <- NULL

# Matrix for a regression-like model
mod <- with(colData(spe_pseudo), model.matrix(~ 0 + age_bin))


# Use pseudoBulkDGE to quickly perform age_bin DE analysis for BayesSpace labels
de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~ 0 + age_bin,
    coef = "age_binInfant",
    row.data = rowData(spe_pseudo)
)


# Save modeling results
saveRDS(de_results, file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
