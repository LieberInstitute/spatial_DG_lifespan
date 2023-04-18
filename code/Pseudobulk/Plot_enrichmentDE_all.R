###############################
# spatial_DG_lifespan project
# Plot DE analysis results
# Anthony Ramnauth, May 18 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

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
    library(ComplexHeatmap)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

################################################################
# Make a .csv list of enriched genes for each BayesSpace Cluster
################################################################

# Make a data frame summary
clust_ML_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_2,
    p_val = modeling_results$enrichment$p_value_2,
    FDR = modeling_results$enrichment$fdr_2
)

clust_ML_enriched_summary <- arrange(clust_ML_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_1 <- file.path(dir_outputs, "cluster_ML_enriched_results")

# Export summary as .csv file
write.csv(clust_ML_enriched_summary,fn_out_1, row.names = FALSE)

# Make a data frame summary
clust_CA3_4_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_4,
    p_val = modeling_results$enrichment$p_value_4,
    FDR = modeling_results$enrichment$fdr_4
)

clust_CA3_4_enriched_summary <- arrange(clust_CA3_4_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_2 <- file.path(dir_outputs, "cluster_CA3&4_enriched_results")

# Export summary as .csv file
write.csv(clust_CA3_4_enriched_summary,fn_out_2, row.names = FALSE)

# Make a data frame summary
clust_SGZ_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_6,
    p_val = modeling_results$enrichment$p_value_6,
    FDR = modeling_results$enrichment$fdr_6
)

clust_SGZ_enriched_summary <- arrange(clust_SGZ_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_3 <- file.path(dir_outputs, "cluster_SGZ_enriched_results")

# Export summary as .csv file
write.csv(clust_SGZ_enriched_summary,fn_out_3, row.names = FALSE)

# Make a data frame summary
clust_GCL_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_7,
    p_val = modeling_results$enrichment$p_value_7,
    FDR = modeling_results$enrichment$fdr_7
)

clust_GCL_enriched_summary <- arrange(clust_GCL_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_4 <- file.path(dir_outputs, "cluster_GCL_enriched_results")

# Export summary as .csv file
write.csv(clust_GCL_enriched_summary,fn_out_4, row.names = FALSE)

################################################
# Plot top enriched genes per BayesSpace Cluster
################################################

# Get mean expression
mat <- assays(spe_pseudo)$logcounts

# filter
#gIndex <- rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
#mat_filter <- mat[gIndex, ] # subset matrix on just those genes.  want to remove lowly expressed genes.

# Extract the p-values
pvals <- modeling_results$enrichment[, 10:18]
rownames(pvals) <- rownames(mat)

# Extract the t-statistics
t_stat <- modeling_results$enrichment[, 1:9]
rownames(t_stat) <- rownames(mat)

# Extract the FDRs
fdrs <- modeling_results$enrichment[, 19:27]
rownames(fdrs) <- rownames(mat)

### pick top 10 genes per cluster:sample
cluster_specific_indices <- mapply(
    function(t, p, f) {
        oo <- order(t, decreasing = TRUE)[1:10]
    },
    as.data.frame(t_stat),
    as.data.frame(pvals),
    as.data.frame(fdrs)
)
cluster_ind <- unique(as.numeric(cluster_specific_indices))

# Add logcounts from indexed from top genes
exprs_heatmap <- assays(spe_pseudo)[[2]][cluster_ind, ]
rownames(exprs_heatmap) <- rowData(spe_pseudo)$gene_name[cluster_ind]
colnames(exprs_heatmap) <- paste("logcount", 1:142, sep = "")

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

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

exprs_heatmap <- scale_rows(exprs_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "enrichment_heatmap_all.pdf"), width = 12, height = 8)
Heatmap(exprs_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace_cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen"),
        BayesSpace_cluster = c("1" = "#5A5156", "2" = "#E4E1E3", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF"))),
    column_title = "Top 10 transcripts from Differential Enrichment of BayesSpace Clusters",
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_order = Bayes_age_order,
    column_split = spe_pseudo$BayesSpace,
    row_split = 10,
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
