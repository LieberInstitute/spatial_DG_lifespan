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

GO_ALL <- compareCluster(ENTREZID~group+age,
    data = clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

# Save ontology objects

save(GO_CC, GO_MF, GO_BP, GO_ALL,
    file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "DG_comp_enrichedGO.Rdata")
)

# directory to save lists
dir_outputs <- here("processed-data", "pseudobulk_spe", "gene_ontologies")
fn_out_1 <- file.path(dir_outputs, "GO_CC")
fn_out_2 <- file.path(dir_outputs, "GO_MF")
fn_out_3 <- file.path(dir_outputs, "GO_BP")

# Export summary as .csv file
write.csv(GO_CC,fn_out_1, row.names = FALSE)
write.csv(GO_MF,fn_out_2, row.names = FALSE)
write.csv(GO_BP,fn_out_3, row.names = FALSE)

#########################################################################################################

# Truncate GO for terms that look interesting

trunc_CC_list <- c(
    "adherens junction",
    "collagen-containing extracellular matrix",
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
    "extracellular matrix organization",
    "'de novo' post-translational protein folding",
    "maintenance of blood-brain barrier",
    "response to ischemia"
)

trunc_ALL_list <- c(
    "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II",
    "synapse pruning",
    "extracellular matrix organization",
    "response to ischemia",
    "lysosomal membrane",
    "gliogenesis",
    "synapse organization",
    "dendrite development",
    "positive regulation of neurogenesis",
    "translational initiation",
    "aerobic respiration",
    "cellular respiration",
    "myeloid cell activation involved in immune response",
    "maintenance of blood-brain barrier",
    "axon development"
)

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_GO.pdf"), width = 10, height = 8)

dotplot(GO_CC, x="group", showCategory = trunc_CC_list) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dotplot(GO_MF, x="group", showCategory = trunc_MF_list) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dotplot(GO_BP, x="group", showCategory = trunc_BP_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dotplot(GO_ALL, x="group", showCategory = trunc_ALL_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()

# Plotting BP separately to adjust sizes & fonts

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_BP_GO.pdf"), width = 12, height = 8)

dotplot(GO_BP, x="group", showCategory = trunc_BP_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_CC_GO.pdf"), width = 10.5, height = 8)

dotplot(GO_CC, x="group", showCategory = trunc_CC_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_MF_GO.pdf"), width = 10.5, height = 8)

dotplot(GO_MF, x="group", showCategory = trunc_MF_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()


pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_ALL_GO.pdf"), width = 10, height = 8)

dotplot(GO_ALL, x="group", showCategory = trunc_ALL_list, label_format = 60) +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey") +
    facet_grid(~factor(age, levels=c('infant', 'teen', 'adult', 'elderly')))

dev.off()

#####################################################################################################################

# Narrow down to the infant GCL since that has the highest # of DEGs

# Load .csv files with the list of DE genes
Infant_GCL_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "InfantvsNonInfant_BayesSpace7_DE.csv"))

infant_GCL <- Infant_GCL_DE_age_results %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_id, gene_name, logFC)

infant_GCL_entrez <- bitr(infant_GCL$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
infant_GCL_entrez <- infant_GCL_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

infant_GCL <- infant_GCL[infant_GCL$gene_id %in% infant_GCL_entrez$ENSEMBL,]

stopifnot(infant_GCL$gene_id == infant_GCL_entrez$ENSEMBL)

infant_GCL$ENTREZID <- infant_GCL_entrez$ENTREZID

infant_GCL$group <- "upregulated"
infant_GCL$group[infant_GCL$logFC < 0] <- "downregulated"
infant_GCL$age <- "infant"

GCL_GO_CC <- compareCluster(ENTREZID~group,
    data = infant_GCL,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

GCL_GO_MF <- compareCluster(ENTREZID~group,
    data = infant_GCL,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

GCL_GO_BP <- compareCluster(ENTREZID~group,
    data = infant_GCL,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

GCL_GO_ALL <- compareCluster(ENTREZID~group,
    data = infant_GCL,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

trunc__GCL_BP_list <- c(
    "neural nucleus development",
    "regulation of ion transmembrane transport",
    "macroautophagy",
    "transition metal ion homeostasis",
    "regulation of transmembrane transporter activity",
    "gliogenesis",
    "regulation of neurogenesis",
    "synapse organization",
    "oligodendrocyte differentiation",
    "vascular endothelial growth factor signaling pathway",
    "myelination",
    "neuron death",
    "extracellular matrix organization",
    "hematopoietic progenitor cell differentiation",
    "dephosphorylation",
    "striated muscle contraction",
    "ossification",
    "regulation of neuron projection development"
)

pdf(file = here::here("plots", "pseudobulked", "infant_GCL_BP_GO.pdf"))

dotplot(GCL_GO_BP, x="group", showCategory = trunc__GCL_BP_list, label_format = 60) +
        theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey")

dev.off()

#################################################################################################################

# GO analysis on infant ML DEGs not found in GCL

# Load .csv files with the list of DE genes
Infant_ML_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "InfantvsNonInfant_BayesSpace2_DE.csv"))

# Narrow down to the ML only DEGs not contained in GCL

Infant_ML_only_gene_list <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "infant_ML_notGCL_DEGs.csv"))

Infant_ML_only_DE_age_results <- Infant_ML_DE_age_results[Infant_ML_DE_age_results$gene_name %in% Infant_ML_only_gene_list$gene_name,]

# format for GO analysis
infant_ML <- Infant_ML_only_DE_age_results %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(gene_id, gene_name, logFC)

infant_ML_entrez <- bitr(infant_ML$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
infant_ML_entrez <- infant_ML_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

infant_ML <- infant_ML[infant_ML$gene_id %in% infant_ML_entrez$ENSEMBL,]

stopifnot(infant_ML$gene_id == infant_ML_entrez$ENSEMBL)

infant_ML$ENTREZID <- infant_ML_entrez$ENTREZID

infant_ML$group <- "upregulated"
infant_ML$group[infant_ML$logFC < 0] <- "downregulated"
infant_ML$age <- "infant"

ML_GO_CC <- compareCluster(ENTREZID~group,
    data = infant_ML,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ML_GO_MF <- compareCluster(ENTREZID~group,
    data = infant_ML,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ML_GO_BP <- compareCluster(ENTREZID~group,
    data = infant_ML,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ML_GO_ALL <- compareCluster(ENTREZID~group,
    data = infant_ML,
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    ont = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

trunc__ML_ALL_list <- c(
    "ATPase-coupled ion transmembrane transporter activity",
    "ATPase activity, coupled to transmembrane movement of ions, rotational mechanism",
    "proton-transporting ATPase activity, rotational mechanism",
    "pyrophosphate hydrolysis−driven proton transmembrane transporter activity",
    "protein refolding",
    "endothelium development",
    "regulation of protein−containing complex assembly",
    "smoothened signaling pathway",
    "branching involved in blood vessel morphogenesis",
    "growth factor binding",
    "structural molecule activity conferring elasticity"
)

pdf(file = here::here("plots", "pseudobulked", "infant_ML_ALL_GO.pdf"), width = 7.5)

dotplot(ML_GO_ALL, x="group", showCategory = trunc__ML_ALL_list, label_format = 60) +
        theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 70, size = 14, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
    scale_color_gradient(low = "black", high = "grey")

dev.off()
