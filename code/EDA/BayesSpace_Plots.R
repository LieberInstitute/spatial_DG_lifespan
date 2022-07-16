###############################
# spatial_DG_lifespan project
# BayesSpace Plots
# Anthony Ramnauth, May 01 2022
###############################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(RColorBrewer)
    library(ggplot2)
    library(Polychrome)
})

# Set up SGE array job to run k=2 to k = 15
# Found in BayesSpaces.sh shell script line -t 2-15
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

# Plot BayesSpace clusters onto tissue
vis_grid_clus(
    spe = spe,
    clustervar = paste0("bayesSpace_harmony_", k),
    pdf = here("plots", "BayesSpace_plots", paste0("vis_grid_clus_BayesSpace_k", k, ".pdf")),
    sort_clust = FALSE,
    colors = setNames(Polychrome::palette36.colors(k), 1:k),
    spatial = FALSE,
    point_size = 2,
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
