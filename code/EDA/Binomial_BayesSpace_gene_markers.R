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
    library(ComplexHeatmap)
    library(circlize)
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
binom_3 <- markers[[3]]
cluster_3_binom_summary <- data.frame(
    gene_name = rownames(binom_3),
    rank = binom_3$Top,
    p_val = binom_3$p.value,
    FDR = binom_3$FDR
)

cluster_3_binom_summary <- cluster_3_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_3 <- file.path(dir_outputs, "cluster_3_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_3_binom_summary,fn_out_3, row.names = FALSE)

# Make a data frame summary
binom_4 <- markers[[4]]
cluster_4_binom_summary <- data.frame(
    gene_name = rownames(binom_4),
    rank = binom_4$Top,
    p_val = binom_4$p.value,
    FDR = binom_4$FDR
)

cluster_4_binom_summary <- cluster_4_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_4 <- file.path(dir_outputs, "cluster_4_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_4_binom_summary,fn_out_4, row.names = FALSE)

# Make a data frame summary
binom_5 <- markers[[5]]
cluster_5_binom_summary <- data.frame(
    gene_name = rownames(binom_5),
    rank = binom_5$Top,
    p_val = binom_5$p.value,
    FDR = binom_5$FDR
)

cluster_5_binom_summary <- cluster_5_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_5 <- file.path(dir_outputs, "cluster_5_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_5_binom_summary,fn_out_5, row.names = FALSE)

# Make a data frame summary
binom_6 <- markers[[6]]
cluster_6_binom_summary <- data.frame(
    gene_name = rownames(binom_6),
    rank = binom_6$Top,
    p_val = binom_6$p.value,
    FDR = binom_6$FDR
)

cluster_6_binom_summary <- cluster_6_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_6 <- file.path(dir_outputs, "cluster_6_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_6_binom_summary,fn_out_6, row.names = FALSE)

# Make a data frame summary
binom_7 <- markers[[7]]
cluster_7_binom_summary <- data.frame(
    gene_name = rownames(binom_7),
    rank = binom_7$Top,
    p_val = binom_7$p.value,
    FDR = binom_7$FDR
)

cluster_7_binom_summary <- cluster_7_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_7 <- file.path(dir_outputs, "cluster_7_binomial_test_results")

# Export summary as .csv file
write.csv(cluster_7_binom_summary,fn_out_7, row.names = FALSE)

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

###########################################################
# Plot log-fold changes for one cluster over other clusters
###########################################################

f1 = colorRamp2(c(-5, 0, 5), c("blue", "#EEEEEE", "red"))

# Selecting cluster 1 since BayesSpace tissue plot looks like ML
interesting_1 <- markers[[1]]
best_set_1 <- interesting_1[interesting_1$Top <= 5, ]
logFCs_1 <- getMarkerEffects(best_set_1)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster1.pdf"))
Heatmap(logFCs_1,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Cluster 1 from binomial test"
    )
dev.off()

# Selecting cluster 4 since BayesSpace tissue plot looks like CA4
interesting_2 <- markers[[2]]
best_set_2 <- interesting_2[interesting_2$Top <= 5, ]
logFCs_2 <- getMarkerEffects(best_set_2)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster2.pdf"))
Heatmap(logFCs_2,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Cluster 2 from binomial test"
    )
dev.off()

# Selecting cluster 4 since BayesSpace tissue plot looks like GCL
interesting_4 <- markers[[4]]
best_set_4 <- interesting_4[interesting_2$Top <= 5, ]
logFCs_4 <- getMarkerEffects(best_set_4)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster4.pdf"))
Heatmap(logFCs_4,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Cluster 4 from binomial test"
    )
dev.off()

# Selecting cluster 8 since BayesSpace tissue plot looks like SGZ
interesting_8 <- markers[[8]]
best_set_8 <- interesting_8[interesting_2$Top <= 5, ]
logFCs_8 <- getMarkerEffects(best_set_8)

pdf(file = here::here("plots", "top_BayesSpace_genes", "logFC_for_Cluster8.pdf"))
Heatmap(logFCs_8,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Cluster 8 from binomial test"
    )
dev.off()

##############
# Add age bins
##############

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

# Limit analysis to only the Dentate Gyrus clusters
# subset spe data based on BayesSpace clusters for DG
spe_DG <- spe[, spe$bayesSpace_harmony_8 %in% c("1", "2", "4", "8")]

###################################
# Test for marker genes of age bins
###################################

age_markers <- findMarkers(spe_DG, groups = spe_DG$age_bin, test = "binom", direction = "up")

# Returns a list with one DataFrame per cluster
age_markers

# directory to save lists
dir_outputs2 <- here("processed-data", "binomial_test_age_bins")
dir.create(dir_plots2, showWarnings = FALSE, recursive = TRUE)

# Make a data frame summary
binom_infant <- age_markers[["Infant"]]
infant_binom_summary <- data.frame(
    gene_name = rownames(binom_infant),
    rank = binom_infant$Top,
    p_val = binom_infant$p.value,
    FDR = binom_infant$FDR
)

infant_binom_summary <- infant_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_inf <- file.path(dir_outputs2, "Infant_binomial_test_results")

# Export summary as .csv file
write.csv(infant_binom_summary,fn_out_inf, row.names = FALSE)

# Make a data frame summary
binom_teen <- age_markers[["Teen"]]
teen_binom_summary <- data.frame(
    gene_name = rownames(binom_teen),
    rank = binom_teen$Top,
    p_val = binom_teen$p.value,
    FDR = binom_teen$FDR
)

teen_binom_summary <- teen_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_teen <- file.path(dir_outputs2, "Teen_binomial_test_results")

# Export summary as .csv file
write.csv(teen_binom_summary,fn_out_teen, row.names = FALSE)

# Make a data frame summary
binom_adult <- age_markers[["Adult"]]
adult_binom_summary <- data.frame(
    gene_name = rownames(binom_adult),
    rank = binom_adult$Top,
    p_val = binom_adult$p.value,
    FDR = binom_adult$FDR
)

adult_binom_summary <- adult_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_adult <- file.path(dir_outputs2, "Adult_binomial_test_results")

# Export summary as .csv file
write.csv(adult_binom_summary,fn_out_adult, row.names = FALSE)

# Make a data frame summary
binom_eld <- age_markers[["Elderly"]]
eld_binom_summary <- data.frame(
    gene_name = rownames(binom_eld),
    rank = binom_eld$Top,
    p_val = binom_eld$p.value,
    FDR = binom_eld$FDR
)

eld_binom_summary <- eld_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_eld <- file.path(dir_outputs2, "Elderly_binomial_test_results")

# Export summary as .csv file
write.csv(eld_binom_summary,fn_out_eld, row.names = FALSE)

###########################################################
# Plot log-fold changes for one age bin over other age bins
###########################################################

# Selecting Infant top gene markers
interesting_inf <- age_markers[["Infant"]]
best_set_inf <- interesting_inf[interesting_inf$Top <= 20, ]
logFCs_inf <- getMarkerEffects(best_set_inf)

pdf(file = here::here("plots", "Binomial_age_bins", "logFC_for_Infant.pdf"))
Heatmap(logFCs_inf,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Infant from binomial test"
    )
dev.off()

# Selecting Teen top gene markers
interesting_te <- age_markers[["Teen"]]
best_set_te <- interesting_te[interesting_te$Top <= 20, ]
logFCs_te <- getMarkerEffects(best_set_te)

pdf(file = here::here("plots", "Binomial_age_bins", "logFC_for_Teen.pdf"))
Heatmap(logFCs_te,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Teen from binomial test"
    )
dev.off()

# Selecting Adult top gene markers
interesting_ad <- age_markers[["Adult"]]
best_set_ad <- interesting_ad[interesting_ad$Top <= 20, ]
logFCs_ad <- getMarkerEffects(best_set_ad)

pdf(file = here::here("plots", "Binomial_age_bins", "logFC_for_Adult.pdf"))
Heatmap(logFCs_ad,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Adult from binomial test"
    )
dev.off()

# Selecting Elderly top gene markers
interesting_el <- age_markers[["Elderly"]]
best_set_el <- interesting_el[interesting_el$Top <= 20, ]
logFCs_el <- getMarkerEffects(best_set_el)

pdf(file = here::here("plots", "Binomial_age_bins", "logFC_for_Elderly.pdf"))
Heatmap(logFCs_el,
    name = "logFC",
    col = f1,
    column_title = "logFC for top genes of Elderly from binomial test"
    )
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
