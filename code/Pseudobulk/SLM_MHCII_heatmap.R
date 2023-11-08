######################################################################################
# spatial_DG_lifespan project
# Plotting genes from parent GO terms of SLM across age groups (Modified for z-scores)
# Anthony Ramnauth, April 26 2023
######################################################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on BayesSpace clusters for DG
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("1", "2", "4", "6", "7")]

bayes_df <- data.frame(spe_pseudo$BayesSpace)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("1", spe_pseudo.BayesSpace) ~ "SLM",
        grepl("2", spe_pseudo.BayesSpace) ~ "ML",
        grepl("4", spe_pseudo.BayesSpace) ~ "CA3&4",
        grepl("6", spe_pseudo.BayesSpace) ~ "SGZ",
        grepl("7", spe_pseudo.BayesSpace) ~ "GCL"
    ))

colData(spe_pseudo)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("SLM", "ML", "CA3&4", "SGZ", "GCL"))

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 27, 18, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 43, 34, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 59, 50, 57, 58, 62, 52, 51, 53, 55, 54,
    80, 76, 77, 79, 72, 65, 75, 66, 73, 74, 78, 68, 67, 69, 71, 70
)

## Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name
# Find relevant GO terms to use

MHCII <- "GO:0019886"
MHCII <- bitr(MHCII, fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
MHCII <- MHCII$SYMBOL

MHCII <- unique(MHCII)

setdiff(MHCII, rownames(spe_pseudo))

MHCII <- MHCII[! MHCII %in% setdiff(MHCII, rownames(spe_pseudo))]

# Add logcounts for all clusters from GO term genes
MHCII_heatmap <- assays(spe_pseudo)[[2]][MHCII, ]
colnames(MHCII_heatmap) <- paste("logcount", 1:80, sep = "")

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "SLM_antigen_MHC_class_II_enrichment_heatmap.pdf"),
    width = 7, height = 6)

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

Heatmap(MHCII_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(spatial_domain = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(spatial_domain = c("SLM" = "black", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    col = col_fun,
    column_title = "antigen processing and presentation of\nexogenous peptide antigen via MHC class II",
    column_order = Bayes_age_order,
    column_split = spe_pseudo$BayesSpace,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = TRUE
    )

dev.off()

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

MHCII_heatmap <- scale_rows(MHCII_heatmap)

pdf(file = here::here("plots", "pseudobulked", "zscores_SLM_antigen_MHC_class_II_enrichment_heatmap.pdf"),
    width = 7, height = 6)

Heatmap(MHCII_heatmap,
    name = "mean\nlog norm counts\n(centered & scaled)",
    top_annotation = HeatmapAnnotation(spatial_domain = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(spatial_domain = c("SLM" = "black", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "antigen processing and presentation of\nexogenous peptide antigen via MHC class II",
    column_order = Bayes_age_order,
    column_split = spe_pseudo$BayesSpace,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = TRUE
    )

dev.off()
