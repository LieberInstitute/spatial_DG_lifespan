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
man_colors <- c("Astro_A" = "yellow3", "Astro_B" = "blue", "Oligo" = "green", "Excit_F" = "black",
    "Tcell" = "magenta", "Mural" = "brown", "Excit_A" = "grey", "OPC_COP" = "green4",
    "Micro" = "cyan2")

pdf(file = here("plots", "Cell_Type_Deconvolution", "Dominant_cell_types_CARD.pdf"), width = 8, height = 6)
vis_clus(
    spe = spe,
    sampleid = "Br1412",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br2706",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br3942",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors

)

vis_clus(
    spe = spe,
    sampleid = "Br5242",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6023",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8195",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8667",
    clustervar = "dominant_cell_type",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8686",
    clustervar = "dominant_cell_type",
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
