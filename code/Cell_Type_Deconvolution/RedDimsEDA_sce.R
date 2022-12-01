##########################################
# spatial_DG_lifespan project
# Explore Reduced Dimensions of sce object
# Anthony Ramnauth, Nov 05 2022
##########################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(scater)
    library(scran)
    library(ggplot2)
    library(PCAtools)
    library(ggnewscale)
    library(spatialLIBD)
    library(sessioninfo)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_harmony.rds"))

# Percentage of variance explained is in the attributes.
percent.var <- attr(reducedDim(sce), "percentVar")

chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow

# Elbow plot of PCs & plot Reduced Dimensions
pdf(file = here::here("sce_plots", "Elbow_plot_sce.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

pdf(file = here::here("sce_plots", "Selected_PCA_plot_sce.pdf"))

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 2),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 3),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 4),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 5),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 6),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 7),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "PCA", ncomponents = c(1, 8),
    colour_by = "Dataset", point_alpha = 0.3
)

dev.off()

pdf(file = here::here("sce_plots", "Selected_Harmony_plot_sce.pdf"))

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 2),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 3),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 4),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 5),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 6),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 7),
    colour_by = "Dataset", point_alpha = 0.3
)

plotReducedDim(sce,
    dimred = "HARMONY", ncomponents = c(1, 8),
    colour_by = "Dataset", point_alpha = 0.3
)

dev.off()




## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
