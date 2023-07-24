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

###############################################################################################
# Make dataframes of ENTREZIDs and logFCs for each age group (ENTREZID seems to work better...)
###############################################################################################
#############################################################################################################
infant <- Infant_DG_DE_age_results %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_id, gene_name, logFC)

infant_entrez <- bitr(infant$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
infant_entrez <- infant_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

infant <- infant[infant$gene_id %in% infant_entrez$ENSEMBL,]

stopifnot(infant$gene_id == infant_entrez$ENSEMBL)

infant$ENTREZID <- infant_entrez$ENTREZID

infant$group <- "upregulated"
infant$group[infant$logFC < 0] <- "downregulated"
infant$age <- "infant"

###########################################################################################################

teen <- Teen_DG_DE_age_results %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_id, gene_name, logFC)

teen_entrez <- bitr(teen$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
teen_entrez <- teen_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

teen <- teen[teen$gene_id %in% teen_entrez$ENSEMBL,]

stopifnot(teen$gene_id == teen_entrez$ENSEMBL)

teen$ENTREZID <- teen_entrez$ENTREZID

teen$group <- "upregulated"
teen$group[teen$logFC < 0] <- "downregulated"
teen$age <- "teen"

##########################################################################################################

adult <- Adult_DG_DE_age_results %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_id, gene_name, logFC)

adult_entrez <- bitr(adult$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
adult_entrez <- adult_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

adult <- adult[adult$gene_id %in% adult_entrez$ENSEMBL,]

stopifnot(adult$gene_id == adult_entrez$ENSEMBL)

adult$ENTREZID <- adult_entrez$ENTREZID

adult$group <- "upregulated"
adult$group[adult$logFC < 0] <- "downregulated"
adult$age <- "adult"

##########################################################################################################

elderly <- Elderly_DG_DE_age_results %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_id, gene_name, logFC)

elderly_entrez <- bitr(elderly$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
elderly_entrez <- elderly_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

elderly <- elderly[elderly$gene_id %in% elderly_entrez$ENSEMBL,]

stopifnot(elderly$gene_id == elderly_entrez$ENSEMBL)

elderly$ENTREZID <- elderly_entrez$ENTREZID

elderly$group <- "upregulated"
elderly$group[elderly$logFC < 0] <- "downregulated"
elderly$age <- "elderly"

##########################################################################################################

clust_compare <- rbind(infant, teen, adult, elderly)

GO_CC <- compareCluster(ENTREZID~group+age,
    data = clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

GO_MF <- compareCluster(ENTREZID~group+age,
    data = clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

GO_BP <- compareCluster(ENTREZID~group+age,
    data = clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

#########################################################################################################

# Truncate GO for terms that look interesting

trunc_CC_list <- c(
    "adherens junction",
    "histone methyltransferase complex",
    "lamellipodium",
    "apical junction complex",
    "NuRD complex",
    "postsynaptic density membrane",
    "methyltransferase complex",
    "lysosomal membrane",
    "transport vesicle",
    "mitochondrial matrix",
    "vacuolar proton−transporting V-type ATPase complex",
    "proteasome complex",
    "MHC class II protein complex",
    "endocytic vesicle",
    "cytoplasmic stress granule",
    "translation preinitiation complex",
    "respirasome",
    "presynapse",
    "stress fiber"
)

trunc_MF_list <- c(
    "transcription coregulator activity",
    "histone binding",
    "histone demethylase activity",
    "ATP-dependent chromatin remodeler activity",
    "miRNA binding",
    "lyase activity",
    "MHC class II protein complex binding",
    "ATPase−coupled ion transmembrane transporter activity",
    "electron transfer activity",
    "translation initiation factor activity",
    "oxidoreduction−driven active transmembrane transporter activity",
    "antioxidant activity",
    "tau-protein kinase activity",
    "immune receptor activity",
    "MHC class II receptor activity"
)

trunc_BP_list <- c(
    "regulation of neuron projection development",
    "regulation of nervous system development",
    "regulation of neurogenesis",
    "synapse organization",
    "histone modification",
    "gliogenesis",
    "neural precursor cell proliferation",
    "axon development",
    "dendrite development",
    "dendritic spine development",
    "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II",
    "aerobic respiration",
    "myeloid cell activation involved in immune response",
    "regulation of translation",
    "cellular respiration",
    "leukocyte mediated immunity",
    "fear response",
    "'de novo' post-translational protein folding",
    "maintenance of blood-brain barrier",
    "response to ischemia"
)

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_GO.pdf"), width = 10, height = 8)

dotplot(GO_CC, x="group", showCategory = trunc_CC_list) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "grey", high = "black") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dotplot(GO_MF, x="group", showCategory = trunc_MF_list) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "grey", high = "black") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dotplot(GO_BP, x="group", showCategory = trunc_BP_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "grey", high = "black") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()

# Plotting BP separately to approximate same size grid as CC & MF

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_BP_GO.pdf"), width = 12, height = 8)

dotplot(GO_BP, x="group", showCategory = trunc_BP_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "grey", high = "black") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_CC_GO.pdf"), width = 10.5, height = 8)

dotplot(GO_CC, x="group", showCategory = trunc_CC_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "grey", high = "black") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()
