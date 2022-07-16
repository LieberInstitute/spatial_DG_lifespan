###############################
# spatial_DG_lifespan project
# Explore Reduced Dimensions
# Anthony Ramnauth, May 02 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(ggplot2)
    library(PCAtools)
    library(ggnewscale)
    library(spatialLIBD)
    library(sessioninfo)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "Dimensions_plots")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Percentage of variance explained is in the attributes.
percent.var <- attr(reducedDim(spe), "percentVar")

chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow

# Elbow plot of PCs & plot Reduced Dimensions
pdf(file = here::here("plots", "Dimensions_plots", "Elbow_plot_spe.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

pdf(file = here::here("plots", "Dimensions_plots", "Selected_PCA_plot_spe.pdf"))
plotReducedDim(spe,
    dimred = "PCA", ncomponents = 4,
    colour_by = "sample_id", point_size = 0.2, point_alpha = 0.5
)
plotReducedDim(spe,
    dimred = "PCA", ncomponents = 4,
    colour_by = "label", point_size = 0.2, point_alpha = 0.5
)
dev.off()

pdf(file = here::here("plots", "Dimensions_plots", "Selected_HARMONY_plot_spe.pdf"))
plotReducedDim(spe,
    dimred = "HARMONY", ncomponents = 4,
    colour_by = "sample_id", point_size = 0.2, point_alpha = 0.5
)
plotReducedDim(spe,
    dimred = "HARMONY", ncomponents = 4,
    colour_by = "label", point_size = 0.2, point_alpha = 0.5
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
