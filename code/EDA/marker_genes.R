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
    library(sessioninfo)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "marker_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Find marker genes
human_markers <-
  c(
    "SNAP25",
    "SLC17A7",
    "GAD1",
    "P2RY12",
    "GAD2",
    "MBP",
    "MOBP",
    "PCP4",
    "RELN",
    "AQP4",
    "CUX2",
    "CCK",
    "HPCAL1",
    "NR4A2",
    "RORB",
    "PROX1",
    "NECAB1",
    "MPPED1",
    "SLC17A6"
  )

# Locate the marker genes
human_markers_search <- rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_name)]

# Plot marker genes on tissue
for (i in human_markers_search) {
  vis_grid_gene(
    spe = spe,
    geneid = i,
    pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
    assayname = "logcounts",
    minCount = -1,
    viridis = TRUE,
    alpha = 0.5,
    sample_order = unique(spe$sample_id),
    point_size = 2
  )
}

# Find neurogenesis marker genes
neurogenesis_markers <-
  c(
    "DCX",
    "SOX11",
    "NEUROD1",
    "SOX2",
    "NEUROG2",
    "MCM2",
    "MKI67",
    "ELAVL2",
    "TACC2",
    "PALLD",
    "NNAT",
    "DPYSL5",
    "DPYSL3",
    "KLF7",
    "BHLHE22"
  )

# Locate the neurogenesis genes
neurogenesis_markers_search <- rowData(spe)$gene_search[match(neurogenesis_markers, rowData(spe)$gene_name)]

# Plot neurogenesis marker genes on tissue
for (i in neurogenesis_markers_search) {
  vis_grid_gene(
    spe = spe,
    geneid = i,
    pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
    assayname = "logcounts",
    minCount = -1,
    viridis = TRUE,
    alpha = 0.5,
    sample_order = unique(spe$sample_id),
    point_size = 2
  )
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
