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

ManualA <- read.csv(file = here("processed-data","spatialLIBD_manual_annotations", "HPC annotations - lex 100622.csv"))

stopifnot(ManualA$spot_name == colnames(spe))

ManualA <- as.vector(ManualA$ManualAnnotation)

spe$ManualAnnotation <- ManualA

spe$ManualAnnotation <- as.factor(spe$ManualAnnotation)

spe$'10x_graphclust' <- as.factor(spe$'10x_graphclust')


pdf(file = here("plots", "QC_plots", "ManualAnnotations.pdf"), width = 8, height = 6)
vis_clus(
    spe = spe,
    sampleid = "Br1412",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br2706",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br3942",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br5242",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br6023",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br8195",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br8667",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br8686",
    clustervar = "ManualAnnotation",
    spatial = FALSE,
    point_size = 2
)

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
