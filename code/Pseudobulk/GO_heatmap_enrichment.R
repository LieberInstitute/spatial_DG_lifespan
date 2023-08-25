#####################################################################################
# spatial_DG_lifespan project
# Plotting genes from parent GO terms of DG across age groups (Modified for z-scores)
# Anthony Ramnauth, April 26 2023
#####################################################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on BayesSpace clusters for DG
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("2", "4", "6", "7")]

bayes_df <- data.frame(spe_pseudo$BayesSpace)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("2", spe_pseudo.BayesSpace) ~ "ML",
        grepl("4", spe_pseudo.BayesSpace) ~ "CA3&4",
        grepl("6", spe_pseudo.BayesSpace) ~ "SGZ",
        grepl("7", spe_pseudo.BayesSpace) ~ "GCL"
    ))

colData(spe_pseudo)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("ML", "CA3&4", "SGZ", "GCL"))

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 27, 18, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 43, 34, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 59, 50, 57, 58, 62, 52, 51, 53, 55, 54
)

## Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Find relevant GO terms to use

cytotrans <- "GO:0019886"
cytotrans <- bitr(cytotrans, fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
cytotrans <- cytotrans$SYMBOL

cytotrans <- unique(cytotrans)

setdiff(cytotrans, rownames(spe_pseudo))

cytotrans <- cytotrans[! cytotrans %in% setdiff(cytotrans, rownames(spe_pseudo))]

# Add logcounts for all clusters from GO term genes
cytotrans_heatmap <- assays(spe_pseudo)[[2]][cytotrans, ]
colnames(cytotrans_heatmap) <- paste("logcount", 1:64, sep = "")

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "antigen_MHC_class_II_enrichment_heatmap.pdf"))

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

Heatmap(cytotrans_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    col = col_fun,
    column_title = "antigen processing and presentation of exogenous peptide antigen via MHC class II",
    column_title_gp = gpar(fontsize = 10),
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE
    )

dev.off()
