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

################################################################
# Make a .csv list of enriched genes for each BayesSpace Cluster
################################################################

# Make a data frame summary
clust_1_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_1,
    p_val = modeling_results$enrichment$p_value_1,
    FDR = modeling_results$enrichment$fdr_1
)

clust_1_enriched_summary <- clust_1_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

clust_1_enriched_summary <- arrange(clust_1_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_1 <- file.path(dir_outputs, "Clust_1_enriched_results")

# Export summary as .csv file
write.csv(clust_1_enriched_summary,fn_out_1, row.names = FALSE)

# Make a data frame summary
clust_2_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_2,
    p_val = modeling_results$enrichment$p_value_2,
    FDR = modeling_results$enrichment$fdr_2
)

clust_2_enriched_summary <- clust_2_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

clust_2_enriched_summary <- arrange(clust_2_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_2 <- file.path(dir_outputs, "Clust_2_enriched_results")

# Export summary as .csv file
write.csv(clust_2_enriched_summary,fn_out_2, row.names = FALSE)

# Make a data frame summary
GCL_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_3,
    p_val = modeling_results$enrichment$p_value_3,
    FDR = modeling_results$enrichment$fdr_3
)

GCL_enriched_summary <- GCL_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

GCL_enriched_summary <- arrange(GCL_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_GCL<- file.path(dir_outputs, "GCL_enriched_results")

# Export summary as .csv file
write.csv(GCL_enriched_summary,fn_out_GCL, row.names = FALSE)

# Make a data frame summary
SGZ_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_4,
    p_val = modeling_results$enrichment$p_value_4,
    FDR = modeling_results$enrichment$fdr_4
)

SGZ_enriched_summary <- SGZ_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SGZ_enriched_summary <- arrange(SGZ_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SGZ<- file.path(dir_outputs, "SGZ_enriched_results")

# Export summary as .csv file
write.csv(SGZ_enriched_summary,fn_out_SGZ, row.names = FALSE)

# Make a data frame summary
CA4_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_5,
    p_val = modeling_results$enrichment$p_value_5,
    FDR = modeling_results$enrichment$fdr_5
)

CA4_enriched_summary <- CA4_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

CA4_enriched_summary <- arrange(CA4_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_CA4<- file.path(dir_outputs, "CA4_enriched_results")

# Export summary as .csv file
write.csv(CA4_enriched_summary,fn_out_CA4, row.names = FALSE)

# Make a data frame summary
CA3_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_6,
    p_val = modeling_results$enrichment$p_value_6,
    FDR = modeling_results$enrichment$fdr_6
)

CA3_enriched_summary <- CA3_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

CA3_enriched_summary <- arrange(CA3_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_CA3<- file.path(dir_outputs, "CA3_enriched_results")

# Export summary as .csv file
write.csv(CA3_enriched_summary,fn_out_CA3, row.names = FALSE)

# Make a data frame summary
ML_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_7,
    p_val = modeling_results$enrichment$p_value_7,
    FDR = modeling_results$enrichment$fdr_7
)

ML_enriched_summary <- ML_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

ML_enriched_summary <- arrange(ML_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_ML<- file.path(dir_outputs, "ML_enriched_results")

# Export summary as .csv file
write.csv(ML_enriched_summary,fn_out_ML, row.names = FALSE)

# Make a data frame summary
clust_8_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_8,
    p_val = modeling_results$enrichment$p_value_8,
    FDR = modeling_results$enrichment$fdr_8
)

clust_8_enriched_summary <- clust_8_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

clust_8_enriched_summary <- arrange(clust_8_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_8<- file.path(dir_outputs, "Clust_8_enriched_results")

# Export summary as .csv file
write.csv(clust_8_enriched_summary,fn_out_8, row.names = FALSE)

################################################
# Plot top enriched genes per BayesSpace Cluster
################################################

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
