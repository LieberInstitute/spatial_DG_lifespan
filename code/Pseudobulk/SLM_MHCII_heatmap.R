######################################################################################
# spatial_DG_lifespan project
# Plotting genes from parent GO terms of SLM across age groups (Modified for z-scores)
# Anthony Ramnauth, April 26 2023
######################################################################################

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

spe_SLM <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "1")]
dim(spe_SLM)

## Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6
)

## Set gene names as row names for easier plotting
rownames(spe_SLM) <- rowData(spe_SLM)$gene_name

# Find relevant GO terms to use

MHCII <- "GO:0019886"
MHCII <- bitr(MHCII, fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
MHCII <- MHCII$SYMBOL

MHCII <- unique(MHCII)

setdiff(MHCII, rownames(spe_SLM))

MHCII <- MHCII[! MHCII %in% setdiff(MHCII, rownames(spe_SLM))]

# Add logcounts for all clusters from GO term genes
MHCII_heatmap <- assays(spe_SLM)[[2]][MHCII, ]
colnames(MHCII_heatmap) <- paste("logcount", 1:16, sep = "")

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

MHCII_heatmap <- scale_rows(MHCII_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "SLM_antigen_MHC_class_II_enrichment_heatmap.pdf"),
    width = 9, height = 6)

Heatmap(MHCII_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = spe_SLM$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "antigen processing and presentation of exogenous peptide antigen via MHC class II",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = TRUE
    )

dev.off()
