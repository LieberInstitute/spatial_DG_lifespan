##################################################
# spatial_DG_lifespan project
# Plots for Dominant Cell Types from cell2location
# Anthony Ramnauth, June 15 2023
##################################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(RColorBrewer)
    library(ggplot2)
    library(Polychrome)
    library(ggspavis)
    library(dplyr)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

cell_props <- as.data.frame(colData(spe)[, c(44:68)])
for (col in 1:ncol(cell_props)){
    colnames(cell_props)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(cell_props)[col])
}

# Find dominant celltype for each spot
p <- cell_props %>%
 rowwise() %>%
 mutate(row_max = names(.)[which.max(c_across(everything()))])

colData(spe)$dominant_cell_types <- p$row_max

# Set Manual colors
cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3",
    "VLMC" = "firebrick1", "Pericyte" = "red2", "Endoth" = "red", "SMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "magenta", "Astro_2" = "yellow3")

pdf(file = here("plots", "Cell_Type_Deconvolution", "Dominant_cell_types_cell2location.pdf"), width = 12, height = 8)

vis_clus(
    spe = spe,
    sampleid = "Br1412",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br2706",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br2720",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br3942",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br5242",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br5699_new",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6023",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6129_new",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6299_new",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6522",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8181",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8195",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8533",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8667",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8686",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8700",
    clustervar = "dominant_cell_types",
    spatial = FALSE,
    point_size = 2,
    colors = cell_colors
)

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
