###########################################################
# spatial_DG_lifespan project
# Volcano plots for senescence markers & BBB markers
# Anthony Ramnauth, June 29 2023
###########################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(sessioninfo)
    library(SingleCellExperiment)
    library(rafalib)
    library(limma)
    library(edgeR)
    library(scran)
    library(EnhancedVolcano)
    library(dplyr)
    library(org.Hs.eg.db)
    library(clusterProfiler)
})

# load modeling results
load(file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_dentategyrus_results.Rdata")
)

# Senescence gene sets

# Get list of gene-set from for senescence markers (from D. Saul, et al., 2022)
Saul_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Saul_2022.csv"))

# Translate from one species to the other using the orthology
SenMayo <- Saul_2022$Gene.human.

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

##############################################
# Volcano Plots of results for elderly age_bin
##############################################

dentate0elderly <- data.frame(
    gene_name = elderly_de_results[[1]]$gene_name,
    logFC = elderly_de_results[[1]]$logFC,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val
)

dentate1elderly <- data.frame(
    gene_name = elderly_de_results[[2]]$gene_name,
    logFC = elderly_de_results[[2]]$logFC,
    adj.P.Val = elderly_de_results[[2]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_eld <- ifelse(
    dentate1elderly$adj.P.Val < 0.05, "#E2C6A7", "#f0e3d6"
)

selected <- c(SenMayo, Casella_2019_up, Casella_2019_down)

## Assigning colors for each groups of highlited genes
keyvals_eld[dentate1elderly$gene_name %in% SenMayo] <- "#789C25"
keyvals_eld[dentate1elderly$gene_name %in% Casella_2019_up] <- "#006164"
keyvals_eld[dentate1elderly$gene_name %in% Casella_2019_down] <- "purple"

## Legend names
names(keyvals_eld)[keyvals_eld == "#E2C6A7"] <- "Adjusted P-value < 0.05"
names(keyvals_eld)[keyvals_eld == "#f0e3d6"] <- "Not significant"
names(keyvals_eld)[keyvals_eld == "#789C25"] <- "SenMayo"
names(keyvals_eld)[keyvals_eld == "#006164"] <- "Casella_2019_up"
names(keyvals_eld)[keyvals_eld == "purple"] <- "Casella_2019_down"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Elderly__senescence_dentategyrus.pdf"))

EnhancedVolcano(dentate1elderly,
    lab = dentate1elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = selected,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    pointSize = c(ifelse(dentate1elderly$gene_name %in% selected, 6, 2)),
    colAlpha = c(ifelse(dentate1elderly$gene_name %in% selected, 1, 0.2)),
    colCustom = keyvals_eld,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Elderly vs. non-Elderly",
    legendPosition = "bottom"
    ) +
    xlim(c(-4, 4)) +
    ylim(c(0, 10))

dev.off()

################################################################################

# BBB maintenance GO term
BBB <- "GO:0035633"
BBB <- bitr(BBB, fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
BBB <- unique(BBB$SYMBOL)

## Colors for the significant and not significant genes
keyvals_eldb <- ifelse(
    dentate1elderly$adj.P.Val < 0.05, "#E2C6A7", "#f0e3d6"
)

keyvals_eldb[dentate1elderly$gene_name %in% BBB] <- "tomato4"

## Legend names
names(keyvals_eldb)[keyvals_eldb == "#E2C6A7"] <- "Adjusted P-value < 0.05"
names(keyvals_eldb)[keyvals_eldb == "#f0e3d6"] <- "Not significant"
names(keyvals_eldb)[keyvals_eldb == "tomato4"] <- "BBB Maintenance"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Elderly__BBB_dentategyrus.pdf"))

EnhancedVolcano(dentate1elderly,
    lab = dentate1elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = BBB,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    pointSize = c(ifelse(dentate1elderly$gene_name %in% BBB, 6, 2)),
    colAlpha = c(ifelse(dentate1elderly$gene_name %in% BBB, 1, 0.2)),
    colCustom = keyvals_eldb,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Elderly vs. non-Elderly",
    legendPosition = "bottom"
    ) +
    xlim(c(-2, 2)) +
    ylim(c(0, 10))

dev.off()

