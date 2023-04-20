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
})

# Create directory for Top Marker Genes plots
dir_outputs <- here::here("processed-data", "top_Spearman_rank_genes_age")
dir.create(dir_outputs, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

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
