#######################################################
# spatial_DG_lifespan project
# GO enrichment & plotting from DE of age groups of GCR
# Anthony Ramnauth, Oct 17 2022
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
})

# Load .csv files with the list of DE genes

Infant_GCR_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "InfantvsNonInfant_GCR_DE.csv"))

Teen_GCR_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "TeenvsNonTeen_GCR_DE.csv"))

Adult_GCR_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "AdultvsNonAdult_GCR_DE.csv"))

Elderly_GCR_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "ElderlyvsNonElderly_NonDentateGyrus_DE.csv"))

# Make dataframes of ENTREZIDs and logFCs for each age group

infant <- Infant_GCR_DE_age_results %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_name, logFC)

infant_entrez <- bitr(infant$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

infant <- infant[infant$gene_name %in% infant_entrez$SYMBOL,]

stopifnot(infant$gene_name == infant_entrez$SYMBOL)

infant <- data.frame(
    ENTREZID = infant_entrez$ENTREZID,
    logFC = infant$logFC
    )

up_infant <- infant[infant$logFC > 0, ]
down_infant <- infant[infant$logFC < 0, ]

teen <- Teen_GCR_DE_age_results %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_name, logFC)

teen_entrez <- bitr(teen$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

teen <- teen[teen$gene_name %in% teen_entrez$SYMBOL,]

stopifnot(teen$gene_name == teen_entrez$SYMBOL)

teen <- data.frame(
    ENTREZID = teen_entrez$ENTREZID,
    logFC = teen$logFC
    )

up_teen <- teen[teen$logFC > 0, ]
down_teen <- teen[teen$logFC < 0, ]

adult <- Adult_GCR_DE_age_results %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_name, logFC)

adult_entrez <- bitr(adult$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

adult <- adult[adult$gene_name %in% adult_entrez$SYMBOL,]

stopifnot(adult$gene_name == adult_entrez$SYMBOL)

adult <- data.frame(
    ENTREZID = adult_entrez$ENTREZID,
    logFC = adult$logFC
    )

up_adult <- adult[adult$logFC > 0, ]
down_adult <- adult[adult$logFC < 0, ]

elderly <- Elderly_GCR_DE_age_results %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_name, logFC)

elderly_entrez <- bitr(elderly$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

elderly <- elderly[elderly$gene_name %in% elderly_entrez$SYMBOL,]

stopifnot(elderly$gene_name == elderly_entrez$SYMBOL)

elderly <- data.frame(
    ENTREZID = elderly_entrez$ENTREZID,
    logFC = elderly$logFC
    )

up_elderly <- elderly[elderly$logFC > 0, ]
down_elderly <- elderly[elderly$logFC < 0, ]

clust_compare <- list(
    up_infant$ENTREZID, down_infant$ENTREZID, up_teen$ENTREZID, down_teen$ENTREZID,
    up_adult$ENTREZID, down_adult$ENTREZID, up_elderly$ENTREZID, down_elderly$ENTREZID
)

names(clust_compare) <- c("Infant_up.reg", "Infant_down.reg", "Teen_up.reg", "Teen_down.reg",
    "Adult_up.reg", "Adult_down.reg", "Elderly_up.reg", "Elderly_down.reg")

# Run compare cluster function

comp_CC <- compareCluster(clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

comp_MF <- compareCluster(clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

comp_BP <- compareCluster(clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE,
)

# Plots for comparing age groups

pdf(file = here::here("plots", "pseudobulked", "age_group_GCR_GO.pdf"), width = 18, height = 14)

dotplot(comp_CC, showCategory = 5, label_format = 60) +
    ggtitle("Top 5 GO Cellular Compartment for Age groups in Granular Cell Region")
dotplot(comp_MF, showCategory = 5, label_format = 60) +
    ggtitle("Top 5 GO Molecular Function for Age groups in Granular Cell Region")
dotplot(comp_BP, showCategory = 5, label_format = 60) +
    ggtitle("Top 5 GO Biological Process for Age groups in Granular Cell Region")

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC,
    showCategory = 5, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 5 GO Cellular Compartment for Age groups in Granular Cell Region")
comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF,
    showCategory = 5, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 5 GO Molecular Function for Age groups in Granular Cell Region")
comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP,
    showCategory = 5, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 5 GO Biological Process for Age groups in Granular Cell Region")

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
