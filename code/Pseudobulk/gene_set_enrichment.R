#############################################
# spatial_DG_lifespan project
# Gene-set enrichment of microglia activation
# Anthony Ramnauth, May 22 2023
#############################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

# Make ortholog list of gene-set from mouse data (Hansruedi Mathys et al., 2017) for microglia states

Mathys_2017 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Mathys_2017.csv"))

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
clust_3_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.3, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_7_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.7, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_6_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.6, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

## Format them appropriately
microglia_geneList <- list(
    early_activated_1 = clust_3_list$ENSEMBL,
    early_activated_2 = clust_7_list$ENSEMBL,
    late_activated = clust_6_list$ENSEMBL
)

enriched_microglia <- gene_set_enrichment(
  microglia_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE
)

gene_set_enrichment_plot(
  enriched_microglia,
  xlabs = unique(enriched_microglia$ID),
  PThresh = 12,
  ORcut = 3,
  enrichOnly = TRUE,
  layerHeights = c(0, seq_len(length(unique(enriched_cluster3$test)))) * 15,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)
