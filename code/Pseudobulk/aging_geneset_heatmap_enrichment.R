#####################################################################################
# spatial_DG_lifespan project
# Plotting genes from parent GO terms of DG across age groups (Modified for z-scores)
# Anthony Ramnauth, April 26 2023
#####################################################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Name the spatial domains for easier plotting
df <-
    data.frame(spe_pseudo$sample_id, spe_pseudo$bayesSpace_harmony_10)

df <- df %>%
    mutate(
        BayesSpace = case_when(
            spe_pseudo.bayesSpace_harmony_10 == 1 ~ "SLM",
            spe_pseudo.bayesSpace_harmony_10 == 2 ~ "ML",
            spe_pseudo.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe_pseudo.bayesSpace_harmony_10 == 5 ~ "SR",
            spe_pseudo.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe_pseudo.bayesSpace_harmony_10 == 7 ~ "GCL",
            spe_pseudo.bayesSpace_harmony_10 == 8 ~ "SL",
            spe_pseudo.bayesSpace_harmony_10 == 9 ~ "CA1",
            spe_pseudo.bayesSpace_harmony_10 == 10 ~ "WM",
        )
    )

colData(spe_pseudo)$BayesSpace <-
    factor(df$BayesSpace, levels = c("SLM", "ML", "CA3&4", "SR", "SGZ",
        "GCL", "SL", "CA1", "WM"))

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 2, 11, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 18, 27, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 34, 43, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 50, 59, 57, 58, 62, 52, 51, 53, 55, 54,
    80, 76, 77, 79, 72, 65, 66, 75, 73, 74, 78, 68, 67, 69, 71, 70,
    96, 92, 93, 95, 88, 81, 82, 91, 89, 90, 94, 84, 83, 85, 87, 86,
    112, 108, 109, 111, 104, 97, 98, 107, 105, 106, 110, 100, 99, 101, 103, 102,
    126, 124, 123, 119, 113, 122, 114, 120, 121, 116, 125, 115, 117, 118,
    142, 138, 139, 141, 134, 127, 128, 137, 135, 136, 140, 130, 129, 131, 133, 132
    )

# Get list of aging gene-set from DEGs of age_bins vs. infant
CAS <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "CAS_list.csv"))

setdiff(CAS$gene_name, rownames(spe_pseudo))

# Add logcounts for all clusters from GO term genes
aging_heatmap <- assays(spe_pseudo)[[2]][CAS$gene_name, ]
colnames(aging_heatmap) <- paste("logcount", 1:142, sep = "")

# Use which(rownames(aging_heatmap) == "gene_name_X") to find which gene names to label in heatmap

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "aging_gene_set_enrichment_heatmap.pdf"))

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

Heatmap(aging_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(spatial_domain = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(spatial_domain = c("SLM" = "#5A5156", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SR" = "#FE00FA",
    "SGZ" = "#1CFFCE", "GCL" = "#B00068", "SL" = "#90AD1C", "CA1" =  "#16FF32", "WM" = "#2ED9FF"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    right_annotation =
        rowAnnotation(foo = anno_mark(at = c(17, 24, 16, 15, 25, 27, 41, 20, 6, 40, 38, 4, 18,
            30, 31, 8, 7, 5, 2, 33, 32, 9, 34, 36, 14, 45, 11, 46),
            labels = c("MAL", "ERMN", "OPALIN", "TMEM144", "EVI2A", "ANLN", "SEMA3B", "CAPN3",
                "CNDP1", "PADI2", "NUPR1", "CNTNAP4", "SLC14A1", "HHATL", "HLA-DRB1", "HLA-DPB1",
                "HLA-DRA", "HLA-DPA1", "CD74", "S100A1", "LAMP5", "QDPR", "TMEM176B", "TMEM176A",
                "HLA-DQB1", "SPP1", "TF", "SEPTIN4", "PPP1R14A"))),
    col = col_fun,
    column_title = "Common Aging Signature genes",
    column_title_gp = gpar(fontsize = 10),
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = FALSE
    )

dev.off()
