##########################################################################
# spatial_DG_lifespan project
# Finding Spearman rank corr for gene expression vs. age for each DG layer
# Anthony Ramnauth, April 20 2023
##########################################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(SpatialExperiment)
    library(BayesSpace)
    library(scran)
    library(scater)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
    library(UpSetR)
    })

# Create directory for Top Marker Genes data
dir_outputs <- here::here("processed-data", "top_Spearman_rank_genes_age")
dir.create(dir_outputs, showWarnings = FALSE, recursive = TRUE)

# Create directory for Top Marker Genes plots
dir_outputs34 <- here::here("plots", "top_Spearman_rank_genes_age")
dir.create(dir_outputs34, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load pseudo-bulked SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Remove Visium spots contaminated with Thalamus using
# gene marker TCF7L2 ENSG00000148737
spe <- spe[, which(logcounts(spe)["ENSG00000148737",] <= 1)]

# order spe observations according to age
spe <- spe[, order(spe$age)]

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

# create separate spe object for each hDG layer

spe_ML <- spe[, which(spe$bayesSpace_harmony_10 == "2")]
dim(spe_ML)

spe_CA3_4 <- spe[, which(spe$bayesSpace_harmony_10 == "4")]
dim(spe_CA3_4)

spe_SGZ <- spe[, which(spe$bayesSpace_harmony_10 == "6")]
dim(spe_SGZ)

spe_GCL <- spe[, which(spe$bayesSpace_harmony_10 == "7")]
dim(spe_GCL)

# Create vector of ages for correlation test

#ML
ages_ML <- seq(0, 100, length = 6974)

#CA3&4
ages_CA3 <- seq(0, 100, length = 10505)

#SGZ
ages_SGZ <- seq(0, 100, length = 5771)

#GCL
ages_GCL <- seq(0, 100, length = 5416)

# calculate the Spearman's rank correlation between each row of the sparse matrix and the ages vector

# ML
#------------------------------------------------------------------------------------------------------------------------------
cor_ML <- apply(logcounts(spe_ML), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML <- sapply(cor_ML, function(x) x$estimate)
p_values_ML <- sapply(cor_ML, function(x) x$p.value)

# combine the results into a data frame
cor__MLdf <- data.frame(gene_id = names(cor_ML), correlation = cor_coeffs_ML, adj.p.value = p_values_ML)

# remove the ".rho" characters from the row names
rownames(cor__MLdf) <- sub("\\.rho$", "", rownames(cor__MLdf))

# order the rows based on the correlation values
cor__MLdf <- cor__MLdf[order(cor__MLdf$correlation, decreasing = TRUE),]

fn_out_1 <- file.path(dir_outputs, "Clust_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__MLdf,fn_out_1, row.names = TRUE)
#------------------------------------------------------------------------------------------------------------------------------
# CA3&4
#------------------------------------------------------------------------------------------------------------------------------
cor_results_CA3_4 <- apply(logcounts(spe_CA3_4), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3 <- sapply(cor_results_CA3_4, function(x) x$estimate)
p_values_CA3 <- sapply(cor_results_CA3_4, function(x) x$p.value)

# combine the results into a data frame
cor_CA3_4df <- data.frame(gene_id = names(cor_results_CA3_4), correlation = cor_coeffs_CA3, adj.p.value = p_values_CA3)

# remove the ".rho" characters from the row names
rownames(cor_CA3_4df) <- sub("\\.rho$", "", rownames(cor_CA3_4df))

# order the rows based on the correlation values
cor_CA3_4df <- cor_CA3_4df[order(cor_CA3_4df$correlation, decreasing = TRUE),]

fn_out_2 <- file.path(dir_outputs, "Clust_CA3&4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor_CA3_4df,fn_out_2, row.names = TRUE)
#------------------------------------------------------------------------------------------------------------------------------
# SGZ
#------------------------------------------------------------------------------------------------------------------------------
cor_results_SGZ <- apply(logcounts(spe_SGZ), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ <- sapply(cor_results_SGZ, function(x) x$estimate)
p_values_SGZ <- sapply(cor_results_SGZ, function(x) x$p.value)

# combine the results into a data frame
cor_SGZdf <- data.frame(gene_id = names(cor_results_SGZ), correlation = cor_coeffs_SGZ, adj.p.value = p_values_SGZ)

# remove the ".rho" characters from the row names
rownames(cor_SGZdf) <- sub("\\.rho$", "", rownames(cor_SGZdf))

# order the rows based on the correlation values
cor_SGZdf <- cor_SGZdf[order(cor_SGZdf$correlation, decreasing = TRUE),]

fn_out_3 <- file.path(dir_outputs, "Clust_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor_SGZdf,fn_out_3, row.names = TRUE)
#------------------------------------------------------------------------------------------------------------------------------
# GCL
#------------------------------------------------------------------------------------------------------------------------------
cor_results_GCL <- apply(logcounts(spe_GCL), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL <- sapply(cor_results_GCL, function(x) x$estimate)
p_values_GCL <- sapply(cor_results_GCL, function(x) x$p.value)

# combine the results into a data frame
cor_GCLdf <- data.frame(gene_id = names(cor_results_GCL), correlation = cor_coeffs_GCL, adj.p.value = p_values_GCL)

# remove the ".rho" characters from the row names
rownames(cor_GCLdf) <- sub("\\.rho$", "", rownames(cor_GCLdf))

# order the rows based on the correlation values
cor_GCLdf <- cor_GCLdf[order(cor_GCLdf$correlation, decreasing = TRUE),]

fn_out_4 <- file.path(dir_outputs, "Clust_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor_GCLdf,fn_out_4, row.names = TRUE)

#################################
# Plot some interesting top genes
#################################

# If starting here load results from Spearman rho and load the pseudo-bulked spe
cor__MLdf <- read.csv(file = here::here("processed-data", "top_Spearman_rank_genes_age",
    "Clust_ML_rho_age_test_results.csv"))

cor_CA3_4df <- read.csv(file = here::here("processed-data", "top_Spearman_rank_genes_age",
    "Clust_CA3&4_rho_age_test_results.csv"))

cor_SGZ_df <- read.csv(file = here::here("processed-data", "top_Spearman_rank_genes_age",
   "Clust_SGZ_rho_age_test_results.csv"))

cor_GCL_df <- read.csv(file = here::here("processed-data", "top_Spearman_rank_genes_age",
    "Clust_GCL_rho_age_test_results.csv"))

# Remove adj.p.value > 0.05
cor__MLdf <- cor__MLdf %>%
    dplyr::filter(adj.p.value < 0.05)

cor_CA3_4df <- cor_CA3_4df %>%
    dplyr::filter(adj.p.value < 0.05)

cor_SGZ_df <- cor_SGZ_df %>%
    dplyr::filter(adj.p.value < 0.05)

cor_GCL_df <- cor_GCL_df %>%
    dplyr::filter(adj.p.value < 0.05)

# Get the head and tail 50 for each list
top_cor_MLdf <- rbind((head(cor__MLdf, 100)), tail(cor__MLdf, 100))
top_cor_CA3_4df <- rbind((head(cor_CA3_4df, 100)), tail(cor_CA3_4df, 100))
top_cor_SGZ_df <- rbind((head(cor_SGZ_df, 100)), tail(cor_SGZ_df, 100))
top_cor_GCL_df <- rbind((head(cor_GCL_df, 100)), tail(cor_GCL_df, 100))

# Set rownames for ease of plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Limit pseudo-bulk data to one hDG layer
pseudo_ML <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "2")]
dim(pseudo_ML)

topML <- top_cor_MLdf$gene_id
topML <- topML[! topML %in% setdiff(topML, rownames(pseudo_ML))]

pseudo_CA3_4 <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "4")]
dim(pseudo_CA3_4)

topCA34 <- top_cor_CA3_4df$gene_id
topCA34 <- topCA34[! topCA34 %in% setdiff(topCA34, rownames(pseudo_CA3_4))]

pseudo_SGZ <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "6")]
dim(pseudo_SGZ)

topSGZ <- top_cor_SGZ_df$gene_id
topSGZ <- topSGZ[! topSGZ %in% setdiff(topSGZ, rownames(pseudo_SGZ))]

pseudo_GCL <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "7")]
dim(pseudo_GCL)

topGCL <- top_cor_GCL_df$gene_id
topGCL <- topGCL[! topGCL %in% setdiff(topGCL, rownames(pseudo_GCL))]

## Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6
)

# Plot

# Add logcounts for all clusters from rho term genes
ML_heatmap <- assays(pseudo_ML)[[2]][topML, ]
colnames(ML_heatmap) <- paste("logcount", 1:16, sep = "")

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

ML_heatmap <- scale_rows(ML_heatmap)

# Add logcounts for all clusters from rho term genes
CA3_4_heatmap <- assays(pseudo_CA3_4)[[2]][topCA34, ]
colnames(CA3_4_heatmap) <- paste("logcount", 1:16, sep = "")

CA3_4_heatmap <- scale_rows(CA3_4_heatmap)

# Add logcounts for all clusters from rho term genes
SGZ_heatmap <- assays(pseudo_SGZ)[[2]][topSGZ, ]
colnames(SGZ_heatmap) <- paste("logcount", 1:16, sep = "")

SGZ_heatmap <- scale_rows(SGZ_heatmap)

# Add logcounts for all clusters from rho term genes
GCL_heatmap <- assays(pseudo_GCL)[[2]][topGCL, ]
colnames(GCL_heatmap) <- paste("logcount", 1:16, sep = "")

GCL_heatmap <- scale_rows(GCL_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "top_Spearman_rank_genes_age", "Spearman_rho_heatmap.pdf"),
    width = 12, height = 22)

Heatmap(ML_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = pseudo_ML$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "ML Top 100 & Bottom 100 Spearman rank correlation heatmap",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = FALSE
    )

Heatmap(CA3_4_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = pseudo_CA3_4$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "CA3&4 Top 100 & Bottom 100 Spearman rank correlation heatmap",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = FALSE
    )

Heatmap(SGZ_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = pseudo_SGZ$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "SGZ Top 100 & Bottom 100 Spearman rank correlation heatmap",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = FALSE
    )

Heatmap(GCL_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = pseudo_GCL$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "GCL Top 100 & Bottom 100 Spearman rank correlation heatmap",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = TRUE,
    cluster_rows = FALSE
    )

dev.off()

# Upset plot of overlapping and unique genes for each layer

pos_corML <- head(cor__MLdf$gene_id, 100)
neg_corML <- tail(cor__MLdf$gene_id, 100)
pos_corCA3_4 <- head(cor_CA3_4df$gene_id, 100)
neg_corCA3_4 <- tail(cor_CA3_4df$gene_id, 100)
pos_cor_SGZ <- head(cor_SGZ_df$gene_id, 100)
neg_cor_SGZ <- tail(cor_SGZ_df$gene_id, 100)
pos_cor_GCL <- head(cor_GCL_df$gene_id, 100)
neg_cor_GCL <- tail(cor_GCL_df$gene_id, 100)

pos_neg <- list(
    pos_rho_ML = pos_corML,
    neg_rho_ML = neg_corML,
    pos_rho_CA3_4 = pos_corCA3_4,
    neg_rho_CA3_4 = neg_corCA3_4,
    pos_rho_SGZ = pos_cor_SGZ,
    neg_rho_SGZ = neg_cor_SGZ,
    pos_rho_GCL = pos_cor_GCL,
    neg_rho_GCL = neg_cor_GCL
)

pdf(file = here::here("plots", "top_Spearman_rank_genes_age", "Spearman_rho_upset.pdf"),width = 15, height = 10)

upset(fromList(pos_neg),
    order.by = "freq",
    sets = c("pos_rho_ML", "neg_rho_ML", "pos_rho_CA3_4", "neg_rho_CA3_4",
        "pos_rho_SGZ", "neg_rho_SGZ", "pos_rho_GCL", "neg_rho_GCL"),
    queries = list(list(query = intersects, params = list("pos_rho_ML",
    "pos_rho_CA3_4", "pos_rho_SGZ", "pos_rho_GCL"), color = "green", active = T),
        list(query = intersects, params = list("neg_rho_ML",
    "neg_rho_CA3_4", "neg_rho_SGZ", "neg_rho_GCL"),
            color = "red", active = T)),
    mainbar.y.label = "# genes",
    sets.x.label = "# genes",
    text.scale = 3
    )

dev.off()

# Extract the unique DG layer genes
pos_genes_ML <- setdiff(pos_neg[[1]], unlist(pos_neg[-1]))
neg_genes_ML <- setdiff(pos_neg[[2]], unlist(pos_neg[-2]))
pos_genes_CA3_4 <- setdiff(pos_neg[[3]], unlist(pos_neg[-3]))
neg_genes_CA3_4 <- setdiff(pos_neg[[4]], unlist(pos_neg[-4]))
pos_genes_SGZ <- setdiff(pos_neg[[5]], unlist(pos_neg[-5]))
neg_genes_SGZ <- setdiff(pos_neg[[6]], unlist(pos_neg[-6]))
pos_genes_GCL <- setdiff(pos_neg[[7]], unlist(pos_neg[-7]))
neg_genes_GCL <- setdiff(pos_neg[[8]], unlist(pos_neg[-8]))

# Create heatmaps for some of these layer unique genes

# Make a pseudo spe for just DG layers

## subset spe data based on BayesSpace clusters for DG
spe_pseudo_DG <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("2", "4", "6", "7")]

bayes_df <- data.frame(spe_pseudo_DG$BayesSpace)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("2", spe_pseudo_DG.BayesSpace) ~ "ML",
        grepl("4", spe_pseudo_DG.BayesSpace) ~ "CA3&4",
        grepl("6", spe_pseudo_DG.BayesSpace) ~ "SGZ",
        grepl("7", spe_pseudo_DG.BayesSpace) ~ "GCL"
    ))

colData(spe_pseudo_DG)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("ML", "CA3&4", "SGZ", "GCL"))

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order_DG <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 27, 18, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 43, 34, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 59, 50, 57, 58, 62, 52, 51, 53, 55, 54
)

## Set gene names as row names for easier plotting
rownames(spe_pseudo_DG) <- rowData(spe_pseudo_DG)$gene_name
#############################################################################
pos_genes_GCL <- pos_genes_GCL[! pos_genes_GCL %in% setdiff(pos_genes_GCL, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
pos_genes_GCL_heatmap <- assays(spe_pseudo_DG)[[2]][pos_genes_GCL, ]
colnames(pos_genes_GCL_heatmap) <- paste("logcount", 1:64, sep = "")

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

pos_genes_GCL_heatmap <- scale_rows(pos_genes_GCL_heatmap)
#############################################################################
neg_genes_GCL <- neg_genes_GCL[! neg_genes_GCL %in% setdiff(neg_genes_GCL, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
neg_genes_GCL_heatmap <- assays(spe_pseudo_DG)[[2]][neg_genes_GCL, ]
colnames(neg_genes_GCL_heatmap) <- paste("logcount", 1:64, sep = "")

neg_genes_GCL_heatmap <- scale_rows(neg_genes_GCL_heatmap)
################################################################################

#############################################################################
pos_genes_ML <- pos_genes_ML[! pos_genes_ML %in% setdiff(pos_genes_ML, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
pos_genes_ML_heatmap <- assays(spe_pseudo_DG)[[2]][pos_genes_ML, ]
colnames(pos_genes_ML_heatmap) <- paste("logcount", 1:64, sep = "")

pos_genes_ML_heatmap <- scale_rows(pos_genes_ML_heatmap)
################################################################################

#############################################################################
neg_genes_ML <- neg_genes_ML[! neg_genes_ML %in% setdiff(neg_genes_ML, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
neg_genes_ML_heatmap <- assays(spe_pseudo_DG)[[2]][neg_genes_ML, ]
colnames(neg_genes_ML_heatmap) <- paste("logcount", 1:64, sep = "")

neg_genes_ML_heatmap <- scale_rows(neg_genes_ML_heatmap)
################################################################################

#############################################################################
pos_genes_CA3_4 <- pos_genes_CA3_4[! pos_genes_CA3_4 %in% setdiff(pos_genes_CA3_4, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
pos_genes_CA3_4_heatmap <- assays(spe_pseudo_DG)[[2]][pos_genes_CA3_4, ]
colnames(pos_genes_CA3_4_heatmap) <- paste("logcount", 1:64, sep = "")

pos_genes_CA3_4_heatmap <- scale_rows(pos_genes_CA3_4_heatmap)
################################################################################

#############################################################################
neg_genes_CA3_4 <- neg_genes_CA3_4[! neg_genes_CA3_4 %in% setdiff(neg_genes_CA3_4, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
neg_genes_CA3_4_heatmap <- assays(spe_pseudo_DG)[[2]][neg_genes_CA3_4, ]
colnames(neg_genes_CA3_4_heatmap) <- paste("logcount", 1:64, sep = "")

neg_genes_CA3_4_heatmap <- scale_rows(neg_genes_CA3_4_heatmap)
################################################################################

#############################################################################
pos_genes_SGZ <- pos_genes_SGZ[! pos_genes_SGZ %in% setdiff(pos_genes_SGZ, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
pos_genes_SGZ_heatmap <- assays(spe_pseudo_DG)[[2]][pos_genes_SGZ, ]
colnames(pos_genes_SGZ_heatmap) <- paste("logcount", 1:64, sep = "")

pos_genes_SGZ_heatmap <- scale_rows(pos_genes_SGZ_heatmap)
################################################################################

#############################################################################
neg_genes_SGZ <- neg_genes_SGZ[! neg_genes_SGZ %in% setdiff(neg_genes_SGZ, rownames(spe_pseudo_DG))]

# Add logcounts for all clusters from GO term genes
neg_genes_SGZ_heatmap <- assays(spe_pseudo_DG)[[2]][neg_genes_SGZ, ]
colnames(neg_genes_SGZ_heatmap) <- paste("logcount", 1:64, sep = "")

neg_genes_SGZ_heatmap <- scale_rows(neg_genes_SGZ_heatmap)
################################################################################

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "top_Spearman_rank_genes_age", "DG_rank_cor_unique_genes_heatmap.pdf"),
    width = 8, height = 10)

Heatmap(pos_genes_GCL_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top GCL unique positve age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(neg_genes_GCL_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top GCL unique negative age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(pos_genes_ML_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top ML unique positive age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(neg_genes_ML_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top ML unique negative age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(pos_genes_CA3_4_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top CA3&4 unique positive age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(neg_genes_CA3_4_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top CA3&4 unique negative age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(pos_genes_SGZ_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top SGZ unique positive age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

Heatmap(neg_genes_SGZ_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo_DG$BayesSpace, age = spe_pseudo_DG$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Top SGZ unique negative age rank correlation genes",
    column_order = Bayes_age_order_DG,
    show_column_names = FALSE,
    column_split = spe_pseudo_DG$BayesSpace,
    show_row_names = TRUE
    )

dev.off()



