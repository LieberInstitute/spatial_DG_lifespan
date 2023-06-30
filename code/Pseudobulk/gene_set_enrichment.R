#############################################
# spatial_DG_lifespan project
# Gene-set enrichment
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

##########################################################################################################

# Activated Microglia gene sets

# Get list of gene-set from mouse data (Hansruedi Mathys et al., 2017) for microglia states

Mathys_2017 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Mathys_2017.csv"))

# Get list of gene-set from human data (Yijing Su et al., 2022) for microglia neuroinflammatory cluster
Su_microg_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Su_microg_2022.csv"))

Su_microg_2022 <- Su_microg_2022[Su_microg_2022$Cluster.ID == "MG1",]

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
clust_3_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.3, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_7_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.7, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_6_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.6, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
MG1 <- bitr(Su_microg_2022$Gene, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

## Format them appropriately
microglia_geneList <- list(
    Mathys_2017_early_activated_1 = clust_3_list$ENSEMBL,
    Mathys_2017_early_activated_2 = clust_7_list$ENSEMBL,
    Mathys_2017_late_activated = clust_6_list$ENSEMBL,
    Su_2022_MG1 = MG1$ENSEMBL
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
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

###############################################################################################################

# Dendritically enriched gene sets

# Get list of gene-set from mouse data (Robert R. Stickels et al., 2020) for dendritically enriched gene sets
Stickels_2020 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Stickels_2020.csv"))

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Translate from one species to the other using the orthology
Dcluster_1 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.1,]
Dcluster_2 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.2,]
Dcluster_3 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.3,]
Dcluster_4 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster4,]

#Get the Ensembl IDs
Dcluster_1 <- bitr(Dcluster_1$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_2 <- bitr(Dcluster_2$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_3 <- bitr(Dcluster_3$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_4 <- bitr(Dcluster_4$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

## Format them appropriately
dendritic_geneList <- list(
    dendritic_enriched_1 = Dcluster_1$ENSEMBL,
    dendritic_enriched_2 = Dcluster_2$ENSEMBL,
    dendritic_enriched_3 = Dcluster_3$ENSEMBL,
    dendritic_enriched_4 = Dcluster_4$ENSEMBL
)

enriched_dendrites <- gene_set_enrichment(
  dendritic_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE
)

gene_set_enrichment_plot(
  enriched_dendrites,
  xlabs = unique(enriched_dendrites$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

###################################################################################

# Senescence gene sets

# Get list of gene-set from for senescence markers (from D. Saul, et al., 2022)
Saul_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Saul_2022.csv"))

# Translate from one species to the other using the orthology
sen <- Saul_2022$Gene.human.

#Get the Ensembl IDs
Sen_cluster <- bitr(sen, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

# Get list of gene-set from for senescence markers (from G. Casella et al., 2019)
Casella_2019_up <- c(
    "TMEM159", "CHPF2", "SLC9A7", "PLOD1", "FAM234B", "DHRS7", "SRPX", "SRPX2", "TNFSF13B", "PDLIM1",
    "ELMOD1", "CCND3", "TMEM30A", "STAT1", "RND3", "TMEM59", "SARAF", "SLCO2B1", "ARRDC4", "PAM",
    "WDR78", "CLSTN2", "WDR63", "NCSTN", "SLC16A14", "GPR155", "CLDN1", "JCAD", "BLCAP", "FILIP1L",
    "TAP1", "TNFRSF10C", "SAMD9L", "SMCO3", "POFUT2", "KIAA1671", "LRP10", "BMS1P9", "MT-TA", "MT-TN",
    "MT-TC", "MT-TY", "DIO2", "MAP4K3-DT", "AC002480.1", "LINC02154", "TM4SF1-AS1", "PTCHD4", "H2AFJ",
    "PURPL"
)

Casella_2019_down <- c(
    "MCUB", "FBL", "HIST1H1D", "HIST1H1A", "FAM129A", "ANP32B", "PARP1", "LBR", "SSRP1", "TMSB15A", "CBS",
    "CDCA7L", "HIST1H1E", "CBX2", "HIST2H2AB", "PTMA", "ITPRIPL1", "AC074135.1", "P16", "P21", "TP53"
)

#Get the Ensembl IDs
Casella_2019_up_cluster <- bitr(Casella_2019_up, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Casella_2019_down_cluster <- bitr(Casella_2019_down, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

## Format them appropriately
senescence_geneList <- list(
    Saul_2022_SenMayo = Sen_cluster,
    Casella_2019_up = Casella_2019_up_cluster
)

senescence_down_geneList <- list(
    Casella_2019_down = Casella_2019_down_cluster
)

enriched_senescence <- gene_set_enrichment(
  senescence_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE
)

gene_set_enrichment_plot(
  enriched_senescence,
  xlabs = unique(enriched_senescence$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

depleted_senescence <- gene_set_enrichment(
  senescence_down_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = TRUE
)

gene_set_enrichment_plot(
  depleted_senescence,
  xlabs = unique(depleted_senescence$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

senescence_all <- rbind(enriched_senescence, depleted_senescence)

gene_set_enrichment_plot(
  senescence_all,
  xlabs = unique(senescence_all$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

######################################################################################
# Get genes from GO term maintenance of blood-brain barrier
BBB <- "GO:0035633"
BBB <- bitr(BBB, fromType="GO", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
BBB <- BBB$ENSEMBL

## Format them appropriately
BBB_geneList <- list(
    BBB_maintenance = BBB
)

enriched_BBB <- gene_set_enrichment(
  BBB_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = TRUE
)

gene_set_enrichment_plot(
  enriched_BBB,
  xlabs = unique(enriched_BBB$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

#######################################################################################

# Reactive Astrocyte gene sets

# Manually input Pan, A1, & A2 gene sets from mouse data (Laura Clarke et al., 2018)

Clarke_2018 <- list(
    PAN = c("Lcn2", "Steap4", "S1pr3", "Timp1", "Hsbp1", "Cxcl10", "Cd44", "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap"),
    A1 = c("C3", "H2-T23", "Serping1", "H2-D1", "Ggta1", "Ligp1", "Gpp2", "Fbln5", "Fkbp5", "Psmb8", "Srgn", "Amigo2"),
    A2 = c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14", "Stat3")
)

# Translate from one species to the other using the orthology
PAN <- orthology[orthology$Column3 %in% Clarke_2018$PAN,]
A1 <- orthology[orthology$Column3 %in% Clarke_2018$A1,]
A2 <- orthology[orthology$Column3 %in% Clarke_2018$A2,]

# Get list of gene-set from human data (Yijing Su et al., 2022) for microglia neuroinflammatory cluster
Su_astro_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Su_astro_2022.csv"))

Su_astro1_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST1",]
Su_astro6_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST6",]
Su_astro_2022 <- rbind(Su_astro1_2022, Su_astro6_2022)

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
PAN_list <- bitr(PAN$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
A1_list <- bitr(A1$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
A2_list <- bitr(A2$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Su_2022_AST1_6 <- bitr(Su_astro_2022$Gene, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

## Format them appropriately
astro_geneList <- list(
    Clarke_2018_PAN = PAN_list$ENSEMBL,
    Clarke_2018_A1 = A1_list$ENSEMBL,
    Clarke_2018_A2 = A2_list$ENSEMBL,
    Su_2022_AST1_AST6 = Su_2022_AST1_6$ENSEMBL
)

enriched_astroglia <- gene_set_enrichment(
  astro_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE
)

gene_set_enrichment_plot(
  enriched_astroglia,
  xlabs = unique(enriched_astroglia$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)

#################################################################################################

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
    Hochgerner_2018_nIPC = nIPC_list$ENSEMBL,
    Hochgerner_2018_NB1 = NB1_list$ENSEMBL,
    Hochgerner_2018_NB2 = NB2_list$ENSEMBL,
    Hao_2022_mouse_macaque_NPC = macaque_list$ENSEMBL,
    Zhou_2022_mouse_human_imGC = imGC_conserved_list$ENSEMBL,
    Zhou_2022_human_imGC = imGC_conserved_list$ENSEMBL
)

enriched_neurogenesis <- gene_set_enrichment(
  neurogenesis_geneList,
  fdr_cut = 0.05,
  modeling_results = modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE
)

gene_set_enrichment_plot(
  enriched_neurogenesis,
  xlabs = unique(enriched_neurogenesis$ID),
  PThresh = 12,
  ORcut = 2,
  enrichOnly = TRUE,
  mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,
    "YlOrRd")))(50)),
  cex = 1.2
)
