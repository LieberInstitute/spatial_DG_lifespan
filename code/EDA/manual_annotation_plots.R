###############################
# spatial_DG_lifespan project
# Plots for Manual Annotations
# Anthony Ramnauth, Oct 07 2022
###############################

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

ManualA <- read.csv(file = here("processed-data","spatialLIBD_manual_annotations",
    "HPC annotations - lex 100622.csv"))

stopifnot(ManualA$spot_name == colnames(spe))

ManualA <- as.vector(ManualA$ManualAnnotation)

spe$ManualAnnotation <- ManualA

spe$ManualAnnotation <- as.factor(spe$ManualAnnotation)

saveRDS(spe, file = here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Spot plots of Manual Annotations

man_colors <- list(ManualAnnotation = Polychrome::palette36.colors(13))
names(man_colors$ManualAnnotation) <- unique(spe$ManualAnnotation)
# Not working so manually assigning names to the colors
man_colors
man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL-CA1" = "#16FF32", "PCL-CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD", "SUB" = "#AA0DFE")

pdf(file = here("plots", "QC_plots", "ManualAnnotations.pdf"), width = 8, height = 6)
vis_clus(
    spe = spe,
    sampleid = "Br1412",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br2706",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br3942",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors

)

vis_clus(
    spe = spe,
    sampleid = "Br5242",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br6023",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8195",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8667",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2,
    colors = man_colors
)

vis_clus(
    spe = spe,
    sampleid = "Br8686",
    clustervar = "ManualAnnotation",
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
