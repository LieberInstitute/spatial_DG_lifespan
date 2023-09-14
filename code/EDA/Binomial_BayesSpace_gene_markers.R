##########################################
# spatial_DG_lifespan project
# Top Marker Genes for BayesSpace Clusters
# Anthony Ramnauth, May 09 2022
##########################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(spatialLIBD)
    library(SpatialExperiment)
    library(BayesSpace)
    library(scran)
    library(scater)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
})

# Create directory for Top Marker Genes plots
dir_plots <- here::here("plots", "top_BayesSpace_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Remove Visium spots contaminated with Thalamus using
# gene marker TCF7L2 ENSG00000148737
spe <- spe[, which(logcounts(spe)["ENSG00000148737",] <= 1)]

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

#################################################
# Test for marker genes across all tissue samples
#################################################

markers <- findMarkers(spe, groups = spe$bayesSpace_harmony_10, test = "binom", direction = "up")

# Returns a list with one DataFrame per cluster
markers

##################################################################
# Make .csv lists of top markers for clusters corresponding to hDG
##################################################################

##-------------------------------------------------------------------

# Make a data frame summary for ML
binom_ML <- markers[[2]]
clust_ML_binom_summary <- data.frame(
    gene_name = rownames(binom_ML),
    rank = binom_ML$Top,
    p_val = binom_ML$p.value,
    FDR = binom_ML$FDR
)

# directory to save lists
dir_outputs <- here("processed-data", "BayesSpace")
fn_out_1 <- file.path(dir_outputs, "Clust_ML_binomial_test_results")

# Export summary as .csv file
write.csv(clust_ML_binom_summary,fn_out_1, row.names = FALSE)

##-------------------------------------------------------------------

# Make a data frame summary for CA3&4
binom_CA3_4 <- markers[[4]]
clust_CA3_4_binom_summary <- data.frame(
    gene_name = rownames(binom_CA3_4),
    rank = binom_CA3_4$Top,
    p_val = binom_CA3_4$p.value,
    FDR = binom_CA3_4$FDR
)

# directory to save lists
fn_out_2 <- file.path(dir_outputs, "Clust_CA3_4_binomial_test_results")

# Export summary as .csv file
write.csv(clust_CA3_4_binom_summary,fn_out_2, row.names = FALSE)

##-------------------------------------------------------------------

# Make a data frame summary for SGZ
binom_SGZ <- markers[[6]]
cluster_SGZ_binom_summary <- data.frame(
    gene_name = rownames(binom_SGZ),
    rank = binom_SGZ$Top,
    p_val = binom_SGZ$p.value,
    FDR = binom_SGZ$FDR
)

# directory to save lists
fn_out_3 <- file.path(dir_outputs, "cluster_SGZ_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_SGZ_binom_summary,fn_out_3, row.names = FALSE)

##-------------------------------------------------------------------

# Make a data frame summary for GCL
binom_GCL <- markers[[7]]
cluster_GCL_binom_summary <- data.frame(
    gene_name = rownames(binom_GCL),
    rank = binom_GCL$Top,
    p_val = binom_GCL$p.value,
    FDR = binom_GCL$FDR
)

# directory to save lists
fn_out_4 <- file.path(dir_outputs, "cluster_GCL_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_GCL_binom_summary,fn_out_4, row.names = FALSE)

##-------------------------------------------------------------------

###############################################################
# Plot log-fold changes for one hDG cluster over other clusters
###############################################################

f1 = colorRamp2(c(-5, 0, 5), c("blue", "#EEEEEE", "red"))

# Selecting cluster 2 since BayesSpace tissue plot looks like ML
interesting_2 <- markers[[2]]
best_set_2 <- interesting_2[interesting_2$Top <= 5, ]
logFCs_2 <- getMarkerEffects(best_set_2)

pdf(file = here::here("plots", "top_BayesSpace_genes", "binomial_logFC_for_ML.pdf"))
Heatmap(logFCs_2,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of ML from binomial test"
    )
dev.off()

# Selecting cluster 4 since BayesSpace tissue plot looks like CA3&4
interesting_4 <- markers[[4]]
best_set_4 <- interesting_4[interesting_4$Top <= 5, ]
logFCs_4 <- getMarkerEffects(best_set_4)

pdf(file = here::here("plots", "top_BayesSpace_genes", "binomial_logFC_for_CA3&4.pdf"))
Heatmap(logFCs_4,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of CA3&4 from binomial test"
    )
dev.off()

# Selecting cluster 6 since BayesSpace tissue plot looks like SGZ
interesting_6 <- markers[[6]]
best_set_6 <- interesting_6[interesting_6$Top <= 5, ]
logFCs_6 <- getMarkerEffects(best_set_6)

pdf(file = here::here("plots", "top_BayesSpace_genes", "binomial_logFC_for_SGZ.pdf"))
Heatmap(logFCs_6,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of SGZ from binomial test"
    )
dev.off()

# Selecting cluster 7 since BayesSpace tissue plot looks like GCL
interesting_7 <- markers[[7]]
best_set_7 <- interesting_7[interesting_7$Top <= 5, ]
logFCs_7 <- getMarkerEffects(best_set_7)

pdf(file = here::here("plots", "top_BayesSpace_genes", "binomial_logFC_for_GCL.pdf"))
Heatmap(logFCs_7,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of GCL from binomial test"
    )
dev.off()
