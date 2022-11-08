#########################################
# spatial_DG_lifespan project
# Plots for Dominant Cell Types from CARD
# Anthony Ramnauth, Oct 18 2022
#########################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(RColorBrewer)
    library(ggplot2)
    library(Polychrome)
    library(ggspavis)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Set Manual colors
man_colors <- c("Oligo_1" = "plum3", "Oligo_2" = "plum4", "Microglia" = "tan2", "Macrocyte" = "tan1",
    "OPC_1" = "goldenrod", "OPC_2" = "goldenrod3", "InN_LAMP5" = "green", "InN_VIP" = "green1",
    "InN_SST" = "green2", "InN_PV" = "green3", "InN_NR2F2" = "green4", "InN_LHX6" = "lawngreen",
    "InN_MEIS2" = "mediumseagreen", "Cajal_Ret" = "black", "Vasc_LM" = "red", "Artl_S_Muscle" = "red1",
    "Pericyte" = "red2", "Endoth" = "red3", "Vasc_S_Muscle" = "red4", "T_cell" = "tan3",
    "Myeloid" = "tan4", "COP" = "goldenrod4", "GC" = "blue", "CA3_N" = "dodgerblue", "EC_N" = "blue1",
    "Mossy" = "blue4", "CA1_N" = "blue2", "SUB_N" = "blue3", "Astro_1" = "magenta", "Astro_2" = "yellow")

pdf(file = here("plots", "Cell_Type_Deconvolution", "Dominant_cell_types_CARD.pdf"), width = 12, height = 8)
vis_clus(
    spe = spe,
    sampleid = "Br1412",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br2706",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br3942",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors

)

vis_clus(
    spe = spe,
    sampleid = "Br5242",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6023",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8195",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8667",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8686",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
