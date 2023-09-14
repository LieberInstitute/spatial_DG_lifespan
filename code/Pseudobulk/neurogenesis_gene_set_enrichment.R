#####################################################
# spatial_DG_lifespan project
# Neurogenesis Gene-set enrichment for each age group
# Anthony Ramnauth, Aug 21 2023
#####################################################

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
})

# Load modeling results
modeling_results_infant <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "infant_modeling_results.rds"))
modeling_results_teen <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "teen_modeling_results.rds"))
modeling_results_adult <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "adult_modeling_results.rds"))
modeling_results_elderly <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "elderly_modeling_results.rds"))


##########################################################################################################

# Neurogenic gene sets

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Get list of neurogenesis gene-set from mouse data (Hochgerner et al., 2018)
Hochgerner_2018 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Hochgerner_2018.csv"))

# Translate from one species to the other using the orthology
Hochgerner_2018_nIPC <- orthology[orthology$Column3 %in% Hochgerner_2018$nIPC,]
Hochgerner_2018_nIPC <- Hochgerner_2018_nIPC$Column1
Hochgerner_2018_NB1 <- orthology[orthology$Column3 %in% Hochgerner_2018$Neuroblast1,]
Hochgerner_2018_NB1 <- Hochgerner_2018_NB1$Column1
Hochgerner_2018_NB2 <- orthology[orthology$Column3 %in% Hochgerner_2018$Neuroblast2,]
Hochgerner_2018_NB2 <- Hochgerner_2018_NB2$Column1

# List taken from Hao et al 2022 for mouse and macaque
Hao_2022 <- c(
    "ID4", "TMEM47", "MLC1", "ALDOC", "SLC1A3", "SLC1A2", "SOX2", "MKI67", "CKS2",
    "CENPF", "SMC4", "TOP2A", "TMPO", "FXYD6", "NNAT", "SOX4", "DPYSL3", "STMN2",
    "CALB2", "SEMA3C"
)

# Get list of imGC gene-set from human-lifespan data (Yi Zhou et al., 2022)
Zhou_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Zhou_2022.csv"))

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
nIPC_list <- bitr(Hochgerner_2018_nIPC, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
NB1_list <- bitr(Hochgerner_2018_NB1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
NB2_list <- bitr(Hochgerner_2018_NB2, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
macaque_list <- bitr(Hao_2022, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
imGC_conserved_list <- bitr(Zhou_2022$Common.genes, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
imGC_human_list <- bitr(Zhou_2022$Human.specific.genes, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

## Format them appropriately
neurogenesis_geneList <- list(
    Mouse_nIPC = nIPC_list$ENSEMBL,
    Mouse_NB1 = NB1_list$ENSEMBL,
    Mouse_NB2 = NB2_list$ENSEMBL,
    Mouse_Macaque_NPC = macaque_list$ENSEMBL,
    Mouse_Human_imGC = imGC_conserved_list$ENSEMBL,
    Human_imGC = imGC_conserved_list$ENSEMBL
)

enriched_infant_neurogenesis <- gene_set_enrichment(
  neurogenesis_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results_infant,
  model_type = names(modeling_results_infant)[2],
  reverse = FALSE
)

enriched_teen_neurogenesis <- gene_set_enrichment(
  neurogenesis_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results_teen,
  model_type = names(modeling_results_teen)[2],
  reverse = FALSE
)

enriched_adult_neurogenesis <- gene_set_enrichment(
  neurogenesis_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results_adult,
  model_type = names(modeling_results_adult)[2],
  reverse = FALSE
)

enriched_elderly_neurogenesis <- gene_set_enrichment(
  neurogenesis_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results_elderly,
  model_type = names(modeling_results_elderly)[2],
  reverse = FALSE
)


# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "gene_set_enrichment", "neurogenesismarkers_GSE.pdf"),
    width = 5, height = 5)

gene_set_enrichment_plot(
  enriched_infant_neurogenesis,
  xlabs = unique(enriched_infant_neurogenesis$ID),
  PThresh = 12,
  ORcut = 3,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

gene_set_enrichment_plot(
  enriched_teen_neurogenesis,
  xlabs = unique(enriched_teen_neurogenesis$ID),
  PThresh = 12,
  ORcut = 3,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

gene_set_enrichment_plot(
  enriched_adult_neurogenesis,
  xlabs = unique(enriched_adult_neurogenesis$ID),
  PThresh = 12,
  ORcut = 3,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

gene_set_enrichment_plot(
  enriched_elderly_neurogenesis,
  xlabs = unique(enriched_elderly_neurogenesis$ID),
  PThresh = 12,
  ORcut = 3,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

dev.off()
