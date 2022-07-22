##########################################
# spatial_DG_lifespan project
# Top Marker Genes for BayesSpace Clusters
# Anthony Ramnauth, May 09 2022
##########################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(SpatialExperiment)
    library(BayesSpace)
    library(scran)
    library(scater)
    library(dplyr)
    library(pheatmap)
})

# Create directory for Top Marker Genes plots
dir_plots <- here::here("plots", "top_BayesSpace_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

spe$bayesSpace_harmony_8 <- as.factor(spe$bayesSpace_harmony_8)

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

#################################################
# Test for marker genes across all tissue samples
#################################################

markers <- findMarkers(spe, groups = spe$bayesSpace_harmony_8, test = "binom", direction = "up")

# Returns a list with one DataFrame per cluster
markers

#################################################
# Make .csv lists of top markers for each cluster
#################################################

# Make a data frame summary
binom_1 <- markers[[1]]
clust_1_binom_summary <- data.frame(
    gene_name = rownames(binom_1),
    rank = binom_1$Top,
    p_val = binom_1$p.value,
    FDR = binom_1$FDR
)

clust_1_binom_summary <- clust_1_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_1 <- file.path(dir_outputs, "Clust_1_binomial_test_results")

# Export summary as .csv file
write.csv(clust_1_binom_summary,fn_out_1, row.names = FALSE)

# Make a data frame summary
binom_2 <- markers[[2]]
clust_2_binom_summary <- data.frame(
    gene_name = rownames(binom_2),
    rank = binom_2$Top,
    p_val = binom_2$p.value,
    FDR = binom_2$FDR
)

clust_2_binom_summary <- clust_2_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_2 <- file.path(dir_outputs, "Clust_2_binomial_test_results")

# Export summary as .csv file
write.csv(clust_2_binom_summary,fn_out_2, row.names = FALSE)

# Make a data frame summary
binom_GCL <- markers[[3]]
GCL_binom_summary <- data.frame(
    gene_name = rownames(binom_GCL),
    rank = binom_GCL$Top,
    p_val = binom_GCL$p.value,
    FDR = binom_GCL$FDR
)

GCL_binom_summary <- GCL_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_GCL <- file.path(dir_outputs, "GCL_binomial_test_results")

# Export summary as .csv file
write.csv(GCL_binom_summary,fn_out_GCL, row.names = FALSE)

# Make a data frame summary
binom_SGZ <- markers[[4]]
SGZ_binom_summary <- data.frame(
    gene_name = rownames(binom_SGZ),
    rank = binom_SGZ$Top,
    p_val = binom_SGZ$p.value,
    FDR = binom_SGZ$FDR
)

SGZ_binom_summary <- SGZ_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_SGZ <- file.path(dir_outputs, "SGZ_binomial_test_results")

# Export summary as .csv file
write.csv(SGZ_binom_summary,fn_out_SGZ, row.names = FALSE)

# Make a data frame summary
binom_CA4 <- markers[[5]]
CA4_binom_summary <- data.frame(
    gene_name = rownames(binom_CA4),
    rank = binom_CA4$Top,
    p_val = binom_CA4$p.value,
    FDR = binom_CA4$FDR
)

CA4_binom_summary <- CA4_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_CA4 <- file.path(dir_outputs, "CA4_binomial_test_results")

# Export summary as .csv file
write.csv(CA4_binom_summary,fn_out_CA4, row.names = FALSE)

# Make a data frame summary
binom_CA3 <- markers[[6]]
CA3_binom_summary <- data.frame(
    gene_name = rownames(binom_CA3),
    rank = binom_CA3$Top,
    p_val = binom_CA3$p.value,
    FDR = binom_CA3$FDR
)

CA3_binom_summary <- CA3_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_CA3 <- file.path(dir_outputs, "CA3_binomial_test_results")

# Export summary as .csv file
write.csv(CA3_binom_summary,fn_out_CA3, row.names = FALSE)

# Make a data frame summary
binom_ML <- markers[[7]]
ML_binom_summary <- data.frame(
    gene_name = rownames(binom_ML),
    rank = binom_ML$Top,
    p_val = binom_ML$p.value,
    FDR = binom_ML$FDR
)

ML_binom_summary <- ML_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_ML <- file.path(dir_outputs, "ML_binomial_test_results")

# Export summary as .csv file
write.csv(ML_binom_summary,fn_out_ML, row.names = FALSE)

# Make a data frame summary
binom_8 <- markers[[8]]
clust_8_binom_summary <- data.frame(
    gene_name = rownames(binom_8),
    rank = binom_8$Top,
    p_val = binom_8$p.value,
    FDR = binom_8$FDR
)

clust_8_binom_summary <- clust_8_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_8 <- file.path(dir_outputs, "Clust_8_binomial_test_results")

# Export summary as .csv file
write.csv(clust_8_binom_summary,fn_out_8, row.names = FALSE)


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 3 since BayesSpace tissue plot looks like GCL
interesting_3 <- markers[[3]]
best_set_3 <- interesting_3[interesting_3$Top <= 5, ]
logFCs_3 <- getMarkerEffects(best_set_3)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster3_GCL.pdf"))
pheatmap(logFCs_3, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-5, 5, length.out = 101),
    main = "logFC for Cluster 3 (GCL) from binomial test")
dev.off()

# Plot log-transformed normalized expression of top genes for GCL cluster
top_genes_3 <- head(rownames(interesting_3))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster3_GCL.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_3, colour_by = "bayesSpace_harmony_8")
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 4 since BayesSpace tissue plot looks like SGZ
interesting_4 <- markers[[4]]
best_set_4 <- interesting_4[interesting_4$Top <= 5, ]
logFCs_4 <- getMarkerEffects(best_set_4)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster4_SGZ.pdf"))
pheatmap(logFCs_4, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 4 (SGZ)")
dev.off()

# Plot log-transformed normalized expression of top genes for SGZ cluster
top_genes_4 <- head(rownames(interesting_4))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster4_SGZ.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_4, colour_by = "bayesSpace_harmony_8")
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 5 since BayesSpace tissue plot looks like CA4
interesting_5 <- markers[[5]]
best_set_5 <- interesting_5[interesting_5$Top <= 5, ]
logFCs_5 <- getMarkerEffects(best_set_5)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster5_CA4.pdf"))
pheatmap(logFCs_5, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 5 (CA4)")
dev.off()

# Plot log-transformed normalized expression of top genes for CA4 cluster
top_genes_5 <- head(rownames(interesting_5))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster5_CA4.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_5, colour_by = "bayesSpace_harmony_8")
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 7 since BayesSpace tissue plot looks like ML
interesting_7 <- markers[[7]]
best_set_7 <- interesting_7[interesting_7$Top <= 5, ]
logFCs_7 <- getMarkerEffects(best_set_7)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster7_ML.pdf"))
pheatmap(logFCs_7, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 7 (ML)")
dev.off()

# Plot log-transformed normalized expression of top genes for ML cluster
top_genes_7 <- head(rownames(interesting_7))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster7_ML.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_7, colour_by = "bayesSpace_harmony_8")
dev.off()

# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 6 since BayesSpace tissue plot looks like CA3
interesting_6 <- markers[[6]]
best_set_6 <- interesting_6[interesting_6$Top <= 5, ]
logFCs_6 <- getMarkerEffects(best_set_6)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster6_CA3.pdf"))
pheatmap(logFCs_6, breaks = seq(-5, 5, length.out = 101), main = "logFC for Cluster 6 (CA3)")
dev.off()

# Plot log-transformed normalized expression of top genes for CA3 cluster
top_genes_6 <- head(rownames(interesting_6))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_cluster6_CA3.pdf"))
plotExpression(spe, x = "bayesSpace_harmony_8", features = top_genes_6, colour_by = "bayesSpace_harmony_8")
dev.off()

#################################################
# Test for marker genes across age bins
#################################################

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        grepl("Br1412", age_df$spe.sample_id) ~ "Teen",
        grepl("Br2706", age_df$spe.sample_id) ~ "Teen",
        grepl("Br3942", age_df$spe.sample_id) ~ "Adult",
        grepl("Br5242", age_df$spe.sample_id) ~ "Elderly",
        grepl("Br6023", age_df$spe.sample_id) ~ "Elderly",
        grepl("Br8195", age_df$spe.sample_id) ~ "Infant",
        grepl("Br8667", age_df$spe.sample_id) ~ "Adult",
        grepl("Br8686", age_df$spe.sample_id) ~ "Infant"
    ))

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

## subset spe data based on age bin
infant_spe <- spe[, spe$age_bin %in% c("Infant")]

############################
# Markers for Infant age bin
############################

infant_markers <- findMarkers(infant_spe,
    groups = infant_spe$bayesSpace_harmony_8,
    test = "binom", direction = "up"
)

# Returns a list with one DataFrame per cluster
infant_markers

# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 3 since BayesSpace tissue plot looks like GCL
infant_interesting_3 <- infant_markers[[3]]
infant_best_set_3 <- infant_interesting_3[infant_interesting_3$Top <= 5, ]
infant_logFCs_3 <- getMarkerEffects(infant_best_set_3)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Infant_Cluster3_GCL.pdf"))
pheatmap(infant_logFCs_3, breaks = seq(-5, 5, length.out = 101), main = "logFC for Infant Cluster 3 (GCL)")
dev.off()

# Plot log-transformed normalized expression of top genes for GCL cluster
infant_top_genes_3 <- head(rownames(infant_interesting_3))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_Infant_cluster3_GCL.pdf"))
plotExpression(infant_spe,
    x = "bayesSpace_harmony_8", features = infant_top_genes_3,
    colour_by = "bayesSpace_harmony_8"
)
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 4 since BayesSpace tissue plot looks like SGZ
infant_interesting_4 <- infant_markers[[4]]
infant_best_set_4 <- infant_interesting_4[infant_interesting_4$Top <= 5, ]
infant_logFCs_4 <- getMarkerEffects(infant_best_set_4)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Infant_Cluster4_SGZ.pdf"))
pheatmap(infant_logFCs_4, breaks = seq(-5, 5, length.out = 101), main = "logFC for Infant Cluster 4 (SGZ)")
dev.off()

# Plot log-transformed normalized expression of top genes for SGZ cluster
infant_top_genes_4 <- head(rownames(infant_interesting_4))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_Infant_cluster4_SGZ.pdf"))
plotExpression(infant_spe,
    x = "bayesSpace_harmony_8", features = infant_top_genes_4,
    colour_by = "bayesSpace_harmony_8"
)
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 5 since BayesSpace tissue plot looks like CA4
infant_interesting_5 <- infant_markers[[5]]
infant_best_set_5 <- infant_interesting_5[infant_interesting_5$Top <= 5, ]
infant_logFCs_5 <- getMarkerEffects(infant_best_set_5)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Infant_Cluster5_CA4.pdf"))
pheatmap(infant_logFCs_5, breaks = seq(-5, 5, length.out = 101), main = "logFC for Infant Cluster 5 (CA4)")
dev.off()

# Plot log-transformed normalized expression of top genes for CA4 cluster
infant_top_genes_5 <- head(rownames(infant_interesting_5))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_Infant_cluster5_CA4.pdf"))
plotExpression(infant_spe,
    x = "bayesSpace_harmony_8", features = infant_top_genes_5,
    colour_by = "bayesSpace_harmony_8"
)
dev.off()


# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 7 since BayesSpace tissue plot looks like ML
infant_interesting_7 <- infant_markers[[7]]
infant_best_set_7 <- infant_interesting_7[infant_interesting_7$Top <= 5, ]
infant_logFCs_7 <- getMarkerEffects(infant_best_set_7)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Infant_Cluster7_ML.pdf"))
pheatmap(infant_logFCs_7, breaks = seq(-5, 5, length.out = 101), main = "logFC for Infant Cluster 7 (ML)")
dev.off()

# Plot log-transformed normalized expression of top genes for ML cluster
infant_top_genes_7 <- head(rownames(infant_interesting_7))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_Infant_cluster7_ML.pdf"))
plotExpression(infant_spe,
    x = "bayesSpace_harmony_8", features = infant_top_genes_7,
    colour_by = "bayesSpace_harmony_8"
)
dev.off()

# Plot log-fold changes for one cluster over all other clusters
# Selecting cluster 6 since BayesSpace tissue plot looks like CA3
infant_interesting_6 <- infant_markers[[6]]
infant_best_set_6 <- infant_interesting_6[infant_interesting_6$Top <= 5, ]
infant_logFCs_6 <- getMarkerEffects(infant_best_set_6)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Infant_Cluster6_CA3.pdf"))
pheatmap(infant_logFCs_6, breaks = seq(-5, 5, length.out = 101), main = "logFC for Infant Cluster 6 (CA3)")
dev.off()

# Plot log-transformed normalized expression of top genes for CA3 cluster
infant_top_genes_6 <- head(rownames(infant_interesting_6))

pdf(file = here::here("plots", "top_BayesSpace_genes", "top_genes_for_Infant_cluster6_CA3.pdf"))
plotExpression(infant_spe,
    x = "bayesSpace_harmony_8", features = infant_top_genes_6,
    colour_by = "bayesSpace_harmony_8"
)
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
