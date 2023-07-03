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
    "RBFOX3", "SNAP25", "SYT1",
    ## Excitatory Neurons
    "SLC17A7",
    # GC
    "PROX1", "CALB1",
    # Mossy Cells
    "ARHGAP24", "DLC1",
    # CA3 PNs
    "CFAP299", "SYN3",
    # CA2 PNs
    "HGF",
    # CA1 PNs
    "ACVR1C", "SYT13",
    # Sub PNs
    "ROBO1", "COL5A2",
    # Cajal-Retzius
    "RELN",
    ## Inhibitory Neurons
    "GAD1", "GAD2",
    # inhibitory subpopulations (Some from Lukas LC & Keri Martinowich 2022-07-22)
    "SST", "KIT", "CALB2", "TAC1", "CNR1", "PVALB", "CORT", "VIP", "NPY",
    "CRHBP", "CCK", "HTR3A", "NR2F2", "LAMP5",
    # Neuropil (Taken from Stickels et al., 2021 & Cajigas et al., 2012, & Muchun Niu et al., 2023)
    "MAP1A", "CAMK2A", "SEMA5A", "SYP",
    # Astrocytes
    "AQP4", "GFAP", "CHRDL1",
    # Oligodendrocytes
    "MOBP",
    # macrophages / microglia
    "CD163", "C3", "PTPRC", "C1QB",
    # OPCs
    "PDGFRA", "VCAN",
    # COP
    "GPR17", "ADAM33",
    # endothelial / mural (RBPMS)
    "CLDN5", "FLT1", "RBPMS", "TTR",
    # T cells
    "SKAP1", "CD247"
)

pdf(file = here::here("plots","BayesSpace_plots", "dotplot_genemarkers.pdf"),
    width = 12, height = 4)

plotDots(spe, group = "bayesSpace_harmony_10", features = features, exprs_values = "logcounts", color = c("blue", "grey", "red"),
    block = "sample_id", center=TRUE, scale=TRUE) +
    scale_y_discrete(limits = features) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) +
    labs(x = "BayesSpace Clusters", y = "Genes", color = "z-score", size = "Proportion of spots") +
    coord_flip()

dev.off()

