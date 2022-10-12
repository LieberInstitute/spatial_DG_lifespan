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
        "GAD2",
        "MBP",
        "MOBP",
        "RELN",
        "AQP4",
        "CCK",
        "HPCAL1",
        "PROX1",
        "NECAB1",
        "MPPED1",
        "SLC17A6",
        "TNNT2",
        "CALB1",
        "GABRQ",
        "THOC3",
        "KCNJ4",
        "SFRP2",
        "APLNR",
        "CD44",
        "TMEM155",
        "SYT13",
        "SLC25A22",
        "SLIT1",
        "SYT4",
        "APOC1",
        "WIF1",
        "MTRNR2L10"
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
        minCount = 0,
        viridis = FALSE,
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
        "BHLHE22",
        "POSTN",
        "BTBD3",
        "CD24",
        "NPDC1",
        "DNM1",
        "NCDN",
        "CAMK2B",
        "NPTX1",
        "RFX3",
        "MALAT1",
        "GFAP",
        "APC"
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
        minCount = 0,
        viridis = FALSE,
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
