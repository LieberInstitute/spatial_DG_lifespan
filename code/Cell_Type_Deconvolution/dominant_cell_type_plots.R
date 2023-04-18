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
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_CARD.rds"))

# Set Manual colors
man_colors <- c("Oligo_1" = "plum", "Oligo_2" = "plum4", "Microglia" = "tan3", "Macrophage" = "tan4",
    "OPC_1" = "goldenrod", "OPC_2" = "goldenrod3", "InN_LAMP5" = "springgreen", "InN_VIP" = "green1",
    "InN_SST" = "springgreen2", "InN_PV" = "green", "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen",
    "InN_MEIS2" = "springgreen3", "InN_NPY" = "green3", "Cajal_Ret" = "black", "Vasc_LM" = "firebrick1",
    "Artl_S_Muscle" = "red1","Pericyte" = "red2", "Endoth" = "red", "Vasc_S_Muscle" = "firebrick",
    "T_cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC_1" = "blue2", "GC_2" = "deepskyblue3",
    "CA3_N" = "deepskyblue", "Mossy" = "blue1", "CA1_d_N" = "blue", "CA1_v_N" = "deepskyblue1",
    "SUB_N" = "blue3", "CA2_N" = "deepskyblue2", "Astro_1" = "yellow", "Astro_2" = "yellow3")

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
