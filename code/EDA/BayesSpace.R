###############################
# spatial_DG_lifespan project
# BayesSpace Spatial Clustering
# Anthony Ramnauth, Apr 29 2022
###############################

# Set up SGE array job to run k=2 to k = 15
# Found in BayesSpaces.sh shell script line -t 2-15

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

# Create directory for BayesSpace plots
dir_plots <- here::here("plots", "BayesSpace_plots")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here::here("processed-data", "BayesSpace_processed_spe", opt$spetype)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(dir_rdata, "clusters_BayesSpace"), showWarnings = FALSE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <-unique(spe$sample_id)
colData(spe)$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
colData(spe)$col <- colData(spe)$array_col

# Run BayesSpace
message("Running spatialCluster()")
Sys.time()
set.seed(12345)
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, platform = "Visium", nrep = 5000)
Sys.time()


spe$bayesSpace_temp<-spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = file.path(dir_rdata, "clusters_BayesSpace")
)

## Visualize BayesSpace results
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))

vis_grid_clus(
    spe = spe,
    clustervar = paste0("BayesSpace_harmony_k", k),
    pdf_file = here("plots", "BayesSpace_plots", paste0("vis_grid_clus_BayesSpace_k",k,".pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 2,
    sample_order = unique(spe$sample_id)
)

saveRDS(spe, file = here::here("processed-data", "BayesSpace_processed_spe", "bayes_spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
