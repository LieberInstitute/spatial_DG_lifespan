###############################
# spatial_DG_lifespan project
# Plotting marker genes
# Anthony Ramnauth, Apr 21 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(ggnewscale)
    library(spatialLIBD)
    library(scater)
    library(scran)
    library(sessioninfo)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "marker_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

## Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

features = c(## Neurons
    "RBFOX3", "SNAP25",
    ## Excitatory Neurons
    "SLC17A7",
    # CA3 PNs
    "KIT", "TRHDE",
    # CA1 PNs
    "MPPED1", "CLMP",
    # GC
    "PROX1", "CALB1",
    ## Inhibitory Neurons
    "GAD1", "GAD2", "SST",
    # Oligodendrocytes
    "MOBP", "MBP",
    # Astrocytes
    "GFAP", "AQP4",
    # Neuropil (Taken from Stickels et al., 2021 & Cajigas et al., 2012, & Muchun Niu et al., 2023)
    "CAMK2A", "MT-ND6",
    # endothelial / mural (RBPMS)
    "CLDN5", "TTR"
    )

pdf(file = here::here("plots","BayesSpace_plots", "dotplot_genemarkers.pdf"),
    width = 6.5, height = 4)

plotDots(spe, group = "bayesSpace_harmony_10", features = features, exprs_values = "logcounts", color = c("white", "red"),
    block = "sample_id") +
    scale_x_discrete(limits = c( "3", "8", "5", "2", "1", "10", "6", "7", "9", "4")) +
    scale_y_discrete(limits = features) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic")) +
    labs(x = "BayesSpace Clusters", y = "Genes", color = "mean\nnorm logcounts", size = "Proportion of spots") +
    coord_flip()

plotDots(spe, group = "bayesSpace_harmony_10", features = features, exprs_values = "logcounts", color = c("blue", "white", "red"),
    block = "sample_id", center = TRUE, scale = TRUE,) +
    scale_x_discrete(limits = c( "3", "8", "5", "2", "1", "10", "6", "7", "9", "4")) +
    scale_y_discrete(limits = features) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic")) +
    labs(x = "BayesSpace Clusters", y = "Genes", color = "mean\nlog norm counts\n(centered & scaled)",
        size = "Proportion of spots") +
    coord_flip()

dev.off()

#######################################################################################################

# Violin plots of markers

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

pdf(file = here::here("plots","BayesSpace_plots", "violinplot_genemarkers.pdf"))

plotExpression(spe, x = "bayesSpace_harmony_10", features = features, colour_by = "bayesSpace_harmony_10") +
    scale_color_manual(values  = bay_colors) +
    theme(plot.title = element_text(face = "italic"))

dev.off()

