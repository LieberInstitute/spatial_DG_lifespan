#########################################################
# spatial_DG_lifespan project
# Plotting top genes from DE of infant DG spatial domains
# Anthony Ramnauth, Oct 11 2023
#########################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
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

# Load .csv that have significant infant DE results & filter for logFCs > 1.5 & < -1.5

infant_vs_non_ML <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "InfantvsNonInfant_BayesSpace2_DE.csv"))

ML_up <- infant_vs_non_ML[infant_vs_non_ML$logFC >= 1.5,]
ML_down <- infant_vs_non_ML[infant_vs_non_ML$logFC <= -1.5,]

infant_vs_non_GCL <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "InfantvsNonInfant_BayesSpace7_DE.csv"))

GCL_up <- infant_vs_non_GCL[infant_vs_non_GCL$logFC >= 1.5,]
GCL_down <- infant_vs_non_GCL[infant_vs_non_GCL$logFC <= -1.5,]

infant_vs_non_SGZ <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "InfantvsNonInfant_BayesSpace6_DE.csv"))

SGZ_up <- infant_vs_non_SGZ[infant_vs_non_SGZ$logFC >= 1.5,]
SGZ_down <- infant_vs_non_SGZ[infant_vs_non_SGZ$logFC <= -1.5,]

infant_vs_non_CA3_4 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "InfantvsNonInfant_BayesSpace4_DE.csv"))

CA3_4_up <- infant_vs_non_CA3_4[infant_vs_non_CA3_4$logFC >= 1.5,]
CA3_4_down <- infant_vs_non_CA3_4[infant_vs_non_CA3_4$logFC <= -1.5,]

# Get top 5 up and down for each

combinedtop <- rbind(ML_up[1:5,], ML_down[1:5,], GCL_up[1:5,], GCL_down[1:5,],
    SGZ_up[1:5,], SGZ_down[1:5,], CA3_4_up[1:5,], CA3_4_down[1:5,])

top_genes <- unique(combinedtop$gene_name)

# Add logcounts for all clusters from GO term genes
top_genes_heatmap <- assays(spe_pseudo)[[2]][top_genes, ]
colnames(top_genes_heatmap) <- paste("logcount", 1:64, sep = "")

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Top_5_DE_infant_vs_noninfant.pdf"), height = 10)

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

Heatmap(top_genes_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    col = col_fun,
    column_title = "Top 5 enriched/depleted infant DE results in DG spatial domains",
    column_title_gp = gpar(fontsize = 10),
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE
    )

dev.off()

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

top_genes_heatmap <- scale_rows(top_genes_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "zscores_Top_5_DE_infant_vs_noninfant.pdf"), height = 10)

Heatmap(top_genes_heatmap,
    name = "mean\nlog norm counts\n(centered & scaled)",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top 5 enriched/depleted infant DE results in DG spatial domains",
    column_title_gp = gpar(fontsize = 10),
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE
    )

dev.off()
