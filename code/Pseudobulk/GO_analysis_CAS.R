#######################################################
# spatial_DG_lifespan project
# GO enrichment & plotting from DE of age_bin vs infant
# Anthony Ramnauth, Sept 6 2023
#######################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(dplyr)
    library(sessioninfo)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(rrvgo)
})

# Load .csv files with the list of DE genes

CAS_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "CAS__BayesSpace_list.csv"))

###############################################################################################
# Make dataframes of ENTREZIDs and logFCs for each age group (ENTREZID seems to work better...)
###############################################################################################
#############################################################################################################

CAS_entrez <- bitr(CAS_results$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

CAS_results <- CAS_results[CAS_results$gene_id %in% CAS_entrez$ENSEMBL,]

stopifnot(CAS_results$gene_id == CAS_entrez$ENSEMBL)

CAS_results$ENTREZID <- CAS_entrez$ENTREZID

CAS_results$group <- "up"
CAS_results$group[CAS_results$sign < 0] <- "down"

####################################################################################################################

GO_ALL <- compareCluster(ENTREZID~group,
    data = CAS_results,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

#########################################################################################################

# Truncate GO for terms that look interesting

trunc_ALL_list <- c(
    "positive regulation of neurogenesis",
    "stem cell differentiation",
    "neuroblast proliferation",
    "Wnt-protein binding",
    "gliogenesis",
    "oligodendrocyte differentiation",
    "Wnt signaling pathway",
    "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II",
    "leukocyte mediated immunity",
    "lymphocyte mediated immunity",
    "transport vesicle membrane",
    "myelin sheath",
    "immune receptor activity"
)

pdf(file = here::here("plots", "pseudobulked", "Select_Aging_GO.pdf"), width = 7, height = 7)

dotplot(GO_ALL, x="group", showCategory = trunc_ALL_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_fill_gradient(low = "black", high = "grey")

dev.off()

