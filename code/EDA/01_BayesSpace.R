################################
# spatial_DG_lifespan project
# BayesSpace Spatial Clustering
# Anthony Ramnauth, Sept 19 2022
################################

# Set up SGE array job to run k=2 to k = 10
# Found in BayesSpaces.sh shell script line -t 2-10

suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(BayesSpace)
    library(ggplot2)
    library(Polychrome)
})

# Create directory for BayesSpace plots
dir_plots <- here::here("plots", "BayesSpace_plots")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <-
    readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

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
        q = k,
        platform = "Visium",
        nrep = 50000
    )
Sys.time()


spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data", "clustering_results"))
