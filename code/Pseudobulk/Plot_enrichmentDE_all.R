###############################
# spatial_DG_lifespan project
# Plot DE analysis results
# Anthony Ramnauth, May 18 2022
###############################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(pheatmap)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data","pseudobulk_spe","modeling_results.rds"))

# Get mean expression
mat <- assays(spe_pseudo)$logcounts

# filter
gIndex = rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter = mat[gIndex, ] #subset matrix on just those genes.  want to remove lowly expressed genes.

# Extract the p-values
pvals <- modeling_results$enrichment[,9:16]
rownames(pvals) = rownames(mat_filter)

# Extract the t-statistics
t_stat <- modeling_results$enrichment[,1:8]
rownames(t_stat) = rownames(mat_filter)

#Extract the FDRs
fdrs <- modeling_results$enrichment[,17:24]
rownames(fdrs) = rownames(mat_filter)

### pick top 10 genes per cluster:sample
cluster_specific_indices = mapply(function(t, p, f) {
  oo = order(t, decreasing = TRUE)[1:10]
},
as.data.frame(t_stat),
as.data.frame(pvals),
as.data.frame(fdrs))
cluster_ind = unique(as.numeric(cluster_specific_indices))

# Add logcounts from indexed from top genes
exprs_heatmap <- assays(spe_pseudo)[[2]][cluster_ind,]
rownames(exprs_heatmap) <- rowData(spe_pseudo)$gene_name[cluster_ind]
colnames(exprs_heatmap) = paste("logcount", 1:64, sep = "")

# Add annotations for pheatmap
cluster_labels <- as.vector(c(rep("Cluster_1", 8), rep("Cluster_2", 8), rep("GCL", 8), rep("SGZ", 8),
    rep("CA4", 8), rep("CA3", 8), rep("ML", 8), rep("Cluster_8", 8)))

annotation_col <- data.frame(BayesSpace = factor(c(cluster_labels)))
rownames(annotation_col) = colnames(exprs_heatmap)
ann_colors = list(BayesSpace = brewer.pal(8, "Set1"))
names(ann_colors$BayesSpace) <- unique(annotation_col$BayesSpace)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots","pseudobulked","enrichment_heatmap_all.pdf"), width = 8, height = 8)
pheatmap(
    exprs_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    color = inferno(20),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    fontsize_row = 9,
    main = "logcounts from enrichment model top 10 genes per cluster",
)
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
