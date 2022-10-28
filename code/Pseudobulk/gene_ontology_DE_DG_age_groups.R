#######################################################
# spatial_DG_lifespan project
# GO enrichment & plotting from DE of age groups of DG
# Anthony Ramnauth, Oct 28 2022
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
    library(pheatmap)
    library(sessioninfo)
})

# Load .csv files with the list of DE genes

Infant_DG_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "InfantvsNonInfant_DentateGyrus_DE.csv"))

Teen_DG_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "TeenvsNonTeen_DentateGyrus_DE.csv"))

Adult_DG_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "AdultvsNonAdult_DentateGyrus_DE.csv"))

Elderly_DG_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "ElderlyvsNonElderly_DentateGyrus_DE.csv"))

############################################################
# Make dataframes of ENTREZIDs and logFCs for each age group
############################################################

infant <- Infant_DG_DE_age_results %>%
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

teen <- Teen_DG_DE_age_results %>%
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

adult <- Adult_DG_DE_age_results %>%
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

elderly <- Elderly_DG_DE_age_results %>%
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

save(comp_CC, comp_MF, comp_BP,
    file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "DG_comp_enrichedGO.Rdata"))

##############################################
# Plot GO comparisons for comparing age groups
##############################################

pdf(file = here::here("plots", "pseudobulked", "Age_group_DG_GO.pdf"), width = 26, height = 16)

dotplot(comp_CC, showCategory = 5, label_format = 90, font.size = 26) +
    ggtitle("Top 5 GO Cellular Compartment for Age groups in Dentate Gyrus") +
    theme(plot.title = element_text(size = 26),
        axis.text.x = element_text(angle = -45))

dotplot(comp_MF, showCategory = 5, label_format = 90, font.size = 26) +
    ggtitle("Top 5 GO Molecular Function for Age groups in Dentate Gyrus") +
    theme(plot.title = element_text(size = 26),
        axis.text.x = element_text(angle = -45))

dotplot(comp_BP, showCategory = 5, label_format = 90, font.size = 26) +
    ggtitle("Top 5 GO Biological Process for Age groups in Dentate Gyrus")+
    theme(plot.title = element_text(size = 26),
        axis.text.x = element_text(angle = -45))

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC,
    showCategory = 5, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 5 GO Cellular Compartment for Age groups in Dentate Gyrus")
comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF,
    showCategory = 5, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 5 GO Molecular Function for Age groups in Dentate Gyrus")
comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP,
    showCategory = 5, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Top 5 GO Biological Process for Age groups in Dentate Gyrus")

dev.off()

######################################################################
# Calculate GO similarities and find parent GO terms for each ontology
######################################################################

# Calculate similarity matrices

CC_vector <- as.vector(comp_CC@compareClusterResult$ID)
MF_vector <- as.vector(comp_MF@compareClusterResult$ID)
BP_vector <- as.vector(comp_BP@compareClusterResult$ID)

CC_simMatrix <- calculateSimMatrix(
    x = CC_vector,
    orgdb = org.Hs.eg.db,
    keytype = "ENTREZID",
    ont = "CC",
    method = "Wang"
)

MF_simMatrix <- calculateSimMatrix(
    x = MF_vector,
    orgdb = org.Hs.eg.db,
    keytype = "ENTREZID",
    ont = "MF",
    method = "Wang"
)

BP_simMatrix <- calculateSimMatrix(
    x = BP_vector,
    orgdb = org.Hs.eg.db,
    keytype = "ENTREZID",
    ont = "BP",
    method = "Wang"
)

# Reduce GO terms

CC_scores <- setNames(-log10(comp_CC@compareClusterResult$qvalue), comp_CC@compareClusterResult$ID)
MF_scores <- setNames(-log10(comp_MF@compareClusterResult$qvalue), comp_MF@compareClusterResult$ID)
BP_scores <- setNames(-log10(comp_BP@compareClusterResult$qvalue), comp_BP@compareClusterResult$ID)

CC_reduced <- reduceSimMatrix(
    simMatrix = CC_simMatrix,
    scores = CC_scores,
    orgdb = org.Hs.eg.db,
    keytype = "ENTREZID"
)

MF_reduced <- reduceSimMatrix(
    simMatrix = MF_simMatrix,
    scores = MF_scores,
    orgdb = org.Hs.eg.db,
    keytype = "ENTREZID"
)

BP_reduced <- reduceSimMatrix(
    simMatrix = BP_simMatrix,
    scores = BP_scores,
    orgdb = org.Hs.eg.db,
    keytype = "ENTREZID"
)

##################################
# Plot SimMatrices & Reduced Terms
##################################

pdf(file = here::here("plots", "pseudobulked", "DG_age_groups_GO_reduction_all.pdf"), width = 20, height = 14)

heatmapPlot(CC_simMatrix,
            CC_reduced,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=12)

heatmapPlot(MF_simMatrix,
            MF_reduced,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=12)

heatmapPlot(BP_simMatrix,
            BP_reduced,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=9)

scatterPlot(CC_simMatrix, CC_reduced)

scatterPlot(MF_simMatrix, MF_reduced)

scatterPlot(BP_simMatrix, BP_reduced)

treemapPlot(CC_reduced)

treemapPlot(MF_reduced)

treemapPlot(BP_reduced)


dev.off()

# Get list of parent GO terms

CC_parent_GO <- data.frame(
    parent_GO = CC_reduced$parent,
    parent_term = CC_reduced$parentTerm
        )
CC_parent_GO <- unique(CC_parent_GO)

MF_parent_GO <- data.frame(
    parent_GO = MF_reduced$parent,
    parent_term = MF_reduced$parentTerm
    )
MF_parent_GO <- unique(MF_parent_GO)

BP_parent_GO <- data.frame(
    parent_GO = BP_reduced$parent,
    parent_term = BP_reduced$parentTerm
        )
BP_parent_GO <- unique(BP_parent_GO)

save(CC_parent_GO, MF_parent_GO, BP_parent_GO,
    file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "DG_comp_parentGOdfs.Rdata"))



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
