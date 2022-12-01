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
cluster_3_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_3,
    p_val = modeling_results$enrichment$p_value_3,
    FDR = modeling_results$enrichment$fdr_3
)

cluster_3_enriched_summary <- cluster_3_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

cluster_3_enriched_summary <- arrange(cluster_3_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_3<- file.path(dir_outputs, "cluster_3_enriched_results")

# Export summary as .csv file
write.csv(cluster_3_enriched_summary,fn_out_3, row.names = FALSE)

# Make a data frame summary
cluster_4_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_4,
    p_val = modeling_results$enrichment$p_value_4,
    FDR = modeling_results$enrichment$fdr_4
)

cluster_4_enriched_summary <- cluster_4_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

cluster_4_enriched_summary <- arrange(cluster_4_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_4<- file.path(dir_outputs, "cluster_4_enriched_results")

# Export summary as .csv file
write.csv(cluster_4_enriched_summary,fn_out_4, row.names = FALSE)

# Make a data frame summary
cluster_5_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_5,
    p_val = modeling_results$enrichment$p_value_5,
    FDR = modeling_results$enrichment$fdr_5
)

cluster_5_enriched_summary <- cluster_5_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

cluster_5_enriched_summary <- arrange(cluster_5_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_5<- file.path(dir_outputs, "cluster_5_enriched_results")

# Export summary as .csv file
write.csv(cluster_5_enriched_summary,fn_out_5, row.names = FALSE)

# Make a data frame summary
cluster_6_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_6,
    p_val = modeling_results$enrichment$p_value_6,
    FDR = modeling_results$enrichment$fdr_6
)

cluster_6_enriched_summary <- cluster_6_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

cluster_6_enriched_summary <- arrange(cluster_6_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_6<- file.path(dir_outputs, "cluster_6_enriched_results")

# Export summary as .csv file
write.csv(cluster_6_enriched_summary,fn_out_6, row.names = FALSE)

# Make a data frame summary
cluster_7_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_7,
    p_val = modeling_results$enrichment$p_value_7,
    FDR = modeling_results$enrichment$fdr_7
)

cluster_7_enriched_summary <- cluster_7_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

cluster_7_enriched_summary <- arrange(cluster_7_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_7<- file.path(dir_outputs, "cluster_7_enriched_results")

# Export summary as .csv file
write.csv(cluster_7_enriched_summary,fn_out_7, row.names = FALSE)

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
gIndex <- rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter <- mat[gIndex, ] # subset matrix on just those genes.  want to remove lowly expressed genes.

# Extract the p-values
pvals <- modeling_results$enrichment[, 9:16]
rownames(pvals) <- rownames(mat_filter)

# Extract the t-statistics
t_stat <- modeling_results$enrichment[, 1:8]
rownames(t_stat) <- rownames(mat_filter)

# Extract the FDRs
fdrs <- modeling_results$enrichment[, 17:24]
rownames(fdrs) <- rownames(mat_filter)

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
colnames(exprs_heatmap) <- paste("logcount", 1:64, sep = "")

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

exprs_heatmap <- scale_rows(exprs_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "enrichment_heatmap_all.pdf"), width = 12, height = 8)
Heatmap(exprs_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace_cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen"),
        BayesSpace_cluster = c("1" = "orangered", "2" = "orange", "3" = "cyan", "4" = "springgreen3",
            "5" = "brown", "6" = "pink", "7" = "yellow", "8" = "slategrey"))),
    column_title = "Top 10 transcripts from Differential Enrichment of BayesSpace Clusters",
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_order = Bayes_age_order,
    column_split = spe_pseudo$BayesSpace,
    row_split = 29,
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
