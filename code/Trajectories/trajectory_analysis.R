###############################
# spatial_DG_lifespan project
# Trajectory analysis
# Anthony Ramnauth, Aug 31 2023
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(slingshot)
    library(scater)
    library(scran)
    library(spatialLIBD)
    library(sessioninfo)
})

# Create directory for  plots
dir_plots <- here::here("plots", "Trajectories")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# order spe observations according to age
spe <- spe[, order(spe$age)]

# Remove Choroid Plexus cluster (screws up variablity estimates)
spe = spe[, which(spe$bayesSpace_harmony_10 != "3")]

spe_sling <- slingshot(spe, cluster = spe$bayesSpace_harmony_10,
    reducedDim = 'PCA', approx_points = 100, omega = TRUE)

pseudo_paths <- slingPseudotime(spe_sling)

head(pseudo_paths)

shared_pseudo <- rowMeans(pseudo_paths, na.rm=TRUE)

colData(spe)$pseudotime <- shared_pseudo

gg <- plotReducedDim(spe_sling, dimred = "UMAP.HARMONY", colour_by = I(shared_pseudo)) +
    labs(x = "UMAP1", y = "UMAP2")

embedded <- embedCurves(spe_sling, "UMAP.HARMONY")
embedded <- slingCurves(embedded)
for (path in embedded) {
    embedded <- data.frame(path$s[path$ord,])
    gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
}

pdf(file = here::here("plots", "Trajectories", "UMAP_trajectories.pdf"))

gg

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$bayesSpace_harmony_10), alpha = 0.01)) +
    geom_point() +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace") +
    theme_bw() +
    theme_classic()

dev.off()

# Plot mean pseudotime on tissue
vis_grid_gene(
    spe = spe,
    geneid = "pseudotime",
    pdf = here::here("plots", "Trajectories", "pseudotime_spotplot.pdf"),
    minCount = 0,
    viridis = TRUE,
    alpha = 0.5,
    point_size = 2,
    spatial = TRUE,
    image_id = "lowres",
    auto_crop = TRUE
    )

# Plot mean pseudotime no background tissue
vis_grid_gene(
    spe = spe,
    geneid = "pseudotime",
    pdf = here::here("plots", "Trajectories", "pseudotime_spotplot_no_tissue.pdf"),
    minCount = 0,
    viridis = TRUE,
    alpha = 1,
    point_size = 2,
    spatial = FALSE,
    auto_crop = TRUE
    )
