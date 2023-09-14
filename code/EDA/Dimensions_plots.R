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
    library(schex)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "Dimensions_plots")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <-
    readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

percent.var <- attr(reducedDim(spe), "percentVar")

chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow

# Elbow plot of PCs & plot Reduced Dimensions
pdf(file = here::here("plots", "Dimensions_plots", "Elbow_plot_spe.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

# Use schex to circumvent overplotting of spots

hex2v1 <- make_hexbin(spe, nbins = 100,
                   dimension_reduction = "PCA", use_dims=c(1,2))

cols <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL_CA1" = "#16FF32", "PCL_CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD",
    "SUB" = "#AA0DFE", "THAL" = "navy")

label_df <- make_hexbin_label(hex2v1, col="bayesSpace_harmony_10")

man_label_df <- make_hexbin_label(hex2v1, col="ManualAnnotation")

scols <- Polychrome::palette36.colors(length(unique(spe$sample_id)))
names(scols) <- sort(unique(spe$sample_id))

label_sdf <- make_hexbin_label(hex2v1, col="sample_id")

pdf(file = here::here("plots", "Dimensions_plots", "PC2vs1_plot_spe.pdf"))

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "sample_id",
    point_size = 0.2,
    point_alpha = 0.5
)

plot_hexbin_meta(hex2v1, col = "sample_id", action = "majority",
                xlab = "PC1", ylab = "PC2", color = scols) +
    labs(fill = "sample_id") +
    theme_bw() +
    theme_classic()

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "bayesSpace_harmony_10",
    point_size = 0.2,
    point_alpha = 0.5
)

plot_hexbin_meta(hex2v1, col = "ManualAnnotation", action = "majority",
                xlab = "PC1", ylab = "PC2", color = man_colors) +
    labs(fill = "Manual\nAnnotations") +
    theme_bw() +
    theme_classic()

plot_hexbin_meta(hex2v1, col = "bayesSpace_harmony_10", action = "majority",
                xlab = "PC1", ylab = "PC2", color = cols) +
    labs(fill = "BayesSpace") +
    theme_bw() +
    theme_classic()

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "sum_gene",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "sum_umi",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "subsets_mito_percent",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "age",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "NBW",
    point_size = 0.2,
    point_alpha = 0.5
)

dev.off()

# Use schex to circumvent overplotting of spots

hex3v1 <- make_hexbin(spe, nbins = 100,
                   dimension_reduction = "PCA", use_dims=c(1,3))

label_df3 <- make_hexbin_label(hex3v1, col="bayesSpace_harmony_10")

pdf(file = here::here("plots", "Dimensions_plots", "PCA3vs1_plot_spe.pdf"))

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "sample_id",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "bayesSpace_harmony_10",
    point_size = 0.2,
    point_alpha = 0.5
)

plot_hexbin_meta(hex3v1, col = "bayesSpace_harmony_10", action = "majority",
                xlab = "PC1", ylab = "PC3", color = cols) +
    labs(fill = "BayesSpace") +
    theme_bw() +
    theme_classic()

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "sum_gene",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "sum_umi",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "subsets_mito_percent",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "age",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 3),
    colour_by = "NBW",
    point_size = 0.2,
    point_alpha = 0.5
)

dev.off()

# Use schex to circumvent overplotting of spots

hex4v1 <- make_hexbin(spe, nbins = 100,
                   dimension_reduction = "PCA", use_dims=c(1,4))

label_df4 <- make_hexbin_label(hex4v1, col="bayesSpace_harmony_10")

pdf(file = here::here("plots", "Dimensions_plots", "PCA4vs1_plot_spe.pdf"))

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "sample_id",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "bayesSpace_harmony_10",
    point_size = 0.2,
    point_alpha = 0.5
)

plot_hexbin_meta(hex4v1, col = "bayesSpace_harmony_10", action = "majority",
                xlab = "PC1", ylab = "PC4", color = cols) +
    labs(fill = "BayesSpace") +
    theme_bw() +
    theme_classic()

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "sum_gene",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "sum_umi",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "subsets_mito_percent",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "age",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = c(1, 4),
    colour_by = "NBW",
    point_size = 0.2,
    point_alpha = 0.5
)

dev.off()

pdf(file = here::here("plots", "Dimensions_plots", "HARMONY2vs1_plot_spe.pdf"))

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "sample_id",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "bayesSpace_harmony_10",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "sum_gene",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "sum_umi",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "subsets_mito_percent",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "age",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = 2,
    colour_by = "NBW",
    point_size = 0.2,
    point_alpha = 0.5
)

dev.off()

pdf(file = here::here("plots", "Dimensions_plots", "HARMONY3vs1_plot_spe.pdf"))

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "sample_id",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "bayesSpace_harmony_10",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "sum_gene",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "sum_umi",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "subsets_mito_percent",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "age",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 3),
    colour_by = "NBW",
    point_size = 0.2,
    point_alpha = 0.5
)

dev.off()

pdf(file = here::here("plots", "Dimensions_plots", "HARMONY4vs1_plot_spe.pdf"))

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "sample_id",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "bayesSpace_harmony_10",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "sum_gene",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "sum_umi",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "subsets_mito_percent",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "age",
    point_size = 0.2,
    point_alpha = 0.5
)

plotReducedDim(
    spe,
    dimred = "HARMONY",
    ncomponents = c(1, 4),
    colour_by = "NBW",
    point_size = 0.2,
    point_alpha = 0.5
)

dev.off()
