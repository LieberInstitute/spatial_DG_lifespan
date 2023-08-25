#############################################################################
# spatial_DG_lifespan project
# Heatmap of genes associated with senescence
# Anthony Ramnauth, June 19 2023
#############################################################################

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
    library(sessioninfo)
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

# Get list of neurogenesis gene-set from mouse data (Hochgerner et al., 2018)
Saul_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Saul_2022.csv"))

sen <- Saul_2022$Gene.human.

sen <- sen[! sen %in%
        setdiff(sen, rownames(spe_pseudo))]

sen_heatmap <- assays(spe_pseudo)[[2]][sen, ]
colnames(sen_heatmap) <- paste("logcount", 1:80, sep = "")


# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Senescence_genemarkers_heatmap.pdf"),
    width = 8, height = 10)

Heatmap(sen_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(spatial_domain = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(spatial_domain = c("SLM" = "black", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Markers for Senescence from SenMayo gene set",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()
