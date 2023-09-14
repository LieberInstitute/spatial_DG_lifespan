############################################################
# spatial_DG_lifespan project
# Cluster Purity check for BayesSpace & Manual Annotations
# Anthony Ramnauth, March 06 2023
############################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(scuttle)
  library(scater)
  library(scran)
  library(dplyr)
  library(PCAtools)
  library(viridis)
  library(bluster)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Check cluster purity for Manual Annotations

pure_spe_man <- neighborPurity(reducedDim(spe, "HARMONY"), spe$ManualAnnotation)
pure_spe_man

pure_data_man <- as.data.frame(pure_spe_man)
pure_data_man$maximum <- factor(pure_data_man$maximum)
pure_data_man$cluster <- spe$ManualAnnotation

pdf(file = here::here("plots","manual_annotations", "Cluster_Purity_ManualAnnotations.pdf"),
    width = 15, height = 8)

ggplot(pure_data_man, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley")

boxplot(split(pure_data_man$purity, pure_data_man$cluster))

dev.off()

# Check cluster purity for BayesSpace k = 10

pure_spe_bayes <- neighborPurity(reducedDim(spe, "HARMONY"), spe$bayesSpace_harmony_10)
pure_spe_bayes

pure_data_bayes <- as.data.frame(pure_spe_bayes)
pure_data_bayes$maximum <- factor(pure_data_bayes$maximum)
pure_data_bayes$cluster <- factor(spe$bayesSpace_harmony_10)

pdf(file = here::here("plots","BayesSpace_plots", "Cluster_Purity_BayesSpace_k10.pdf"),
    width = 15, height = 8)

ggplot(pure_data_bayes, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley")

boxplot(split(pure_data_bayes$purity, pure_data_bayes$cluster))

dev.off()
