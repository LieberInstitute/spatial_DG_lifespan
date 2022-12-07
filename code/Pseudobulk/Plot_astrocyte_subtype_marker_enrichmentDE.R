###########################################################
# spatial_DG_lifespan project
# Plot DEG for Astrocyte subtype entrichment using markers
# from Karpf, J., et al. (2022) Nat. Neuro.
# Anthony Ramnauth, Dec 07 2022
###########################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Vector of GO "Regulation of Neurogenesis" term genes
astro_genes <- as.character(c(
    "GFAP", "SLC1A3", "SLC1A2", "GLUL", "GRM3", "TNC", "ID3"
))

# Add logcounts for all clusters from GO term genes
astro_heatmap <- assays(spe_pseudo)[[2]][astro_genes, ]
colnames(astro_heatmap) <- paste("logcount", 1:64, sep = "")

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
6, 8, 1, 2, 7, 3, 4, 5,
14, 16, 9, 10, 15, 11, 12, 13,
22, 24, 17, 18, 23, 19, 20, 21,
30, 32, 25, 26, 31, 27, 28, 29,
38, 40, 33, 34, 39, 35, 36, 37,
46, 48, 41, 42, 47, 43, 44, 45,
54, 56, 49, 50, 55, 51, 52, 53,
62, 64, 57, 58, 63, 59, 60, 61
    )

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

astro_heatmap <- scale_rows(astro_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "astrocyte_subtypes_enrichment_heatmap_all.pdf"),
    width = 12, height = 8)
Heatmap(astro_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace_cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen"),
        BayesSpace_cluster = c("1" = "orangered", "2" = "orange", "3" = "cyan", "4" = "springgreen3",
            "5" = "brown", "6" = "pink", "7" = "yellow", "8" = "slategrey"))),
    column_title = "Astrocyte subtype gene marker enrichment within BayesSpace clusters",
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_order = Bayes_age_order,
    column_split = spe_pseudo$BayesSpace,
    row_split = 2,
    row_title = NULL,
    row_names_gp = gpar(fontsize = 7)
    )
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
