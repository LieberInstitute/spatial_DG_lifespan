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
    "NECAB1", "KIT", "CCK",
    # GC
    "PROX1", "CALB1",
    # CA1 PNs
    "MPPED1", "CLMP",
    # Sub PNs
    "COL5A2", "ROBO1",
    ## Inhibitory Neurons
    "GAD1", "GAD2",
    # inhibitory subpopulations (Some from Lukas LC & Keri Martinowich 2022-07-22)
    "SST",
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
    width = 12, height = 4)

plotDots(spe, group = "bayesSpace_harmony_10", features = features, exprs_values = "logcounts", color = c("blue", "grey", "red"),
    block = "sample_id", center=TRUE, scale=TRUE) +
    scale_x_discrete(limits = c( "3", "8", "5", "2", "1", "10", "6", "9", "7", "4")) +
    scale_y_discrete(limits = features) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1, face = "italic")) +
    labs(x = "BayesSpace Clusters", y = "Genes", color = "z-score", size = "Proportion of spots") +
    coord_flip()

dev.off()

