#######################################################################
# spatial_DG_lifespan project
# Spearman rank corr for senescence expression vs. age for SLM WM, & DG
# Anthony Ramnauth, May 22 2023
#######################################################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(SpatialExperiment)
    library(BayesSpace)
    library(dplyr)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(ComplexHeatmap)
    })

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Create directory for Top Marker Genes data
dir_outputs <- here::here("processed-data", "top_Spearman_rank_genes_age")
dir.create(dir_outputs, showWarnings = FALSE, recursive = TRUE)

# Remove Visium spots contaminated with Thalamus using
# gene marker TCF7L2 ENSG00000148737
spe <- spe[, which(logcounts(spe)["ENSG00000148737",] <= 1)]

# order spe observations according to age
spe <- spe[, order(spe$age)]

# create separate spe object for each BayesSpace spatial domain of interest
spe_SLM <- spe[, which(spe$bayesSpace_harmony_10 == "1")]
dim(spe_SLM)

spe_ML <- spe[, which(spe$bayesSpace_harmony_10 == "2")]
dim(spe_ML)

spe_CA3_4 <- spe[, which(spe$bayesSpace_harmony_10 == "4")]
dim(spe_CA3_4)

spe_SGZ <- spe[, which(spe$bayesSpace_harmony_10 == "6")]
dim(spe_SGZ)

spe_GCL <- spe[, which(spe$bayesSpace_harmony_10 == "7")]
dim(spe_GCL)

# Get list of gene-set from mouse data (D. Saul et al., 2022) for microglia states
Saul_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Saul_2022.csv"))

sen <- Saul_2022$Gene.human.

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
sen <- bitr(sen, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
sen <- sen$ENSEMBL

#SLM
ages_SLM <- spe_SLM$age

#ML
ages_ML <- spe_ML$age

#CA3&4
ages_CA3 <- spe_CA3_4$age

#SGZ
ages_SGZ <- spe_SGZ$age

#GCL
ages_GCL <- spe_GCL$age

# Subset for Senescence clusters
spe_SLM_sen <- spe_SLM[which((rownames(spe_SLM)) %in% sen)]
# Set gene names as row names for easier plotting
rownames(spe_SLM_sen) <- rowData(spe_SLM_sen)$gene_name

spe_ML_sen <- spe_ML[which((rownames(spe_ML)) %in% sen)]
# Set gene names as row names for easier plotting
rownames(spe_ML_sen) <- rowData(spe_ML_sen)$gene_name

spe_CA3_sen <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% sen)]
# Set gene names as row names for easier plotting
rownames(spe_CA3_sen) <- rowData(spe_CA3_sen)$gene_name

spe_SGZ_sen <- spe_SGZ[which((rownames(spe_SGZ)) %in% sen)]
# Set gene names as row names for easier plotting
rownames(spe_SGZ_sen) <- rowData(spe_SGZ_sen)$gene_name

spe_GCL_sen <- spe_GCL[which((rownames(spe_GCL)) %in% sen)]
# Set gene names as row names for easier plotting
rownames(spe_GCL_sen) <- rowData(spe_GCL_sen)$gene_name

#Calculat Spearman rank corr
cor_SLM_sen <- apply(logcounts(spe_SLM_sen), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_sen <- sapply(cor_SLM_sen, function(x) x$estimate)
p_values_SLM_sen <- sapply(cor_SLM_sen, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_sendf <- data.frame(gene_id = names(cor_SLM_sen), correlation = cor_coeffs_SLM_sen, adj.p.value = p_values_SLM_sen)

# remove the ".rho" characters from the row names
rownames(cor__SLM_sendf) <- sub("\\.rho$", "", rownames(cor__SLM_sendf))

hist(cor__SLM_sendf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Senescence gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_sendf <- cor__SLM_sendf[order(cor__SLM_sendf$correlation, decreasing = TRUE),]

fn_out_1 <- file.path(dir_outputs, "Senescence_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_sendf,fn_out_1, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Calculat Spearman rank corr
cor_ML_sen <- apply(logcounts(spe_ML_sen), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_sen <- sapply(cor_ML_sen, function(x) x$estimate)
p_values_ML_sen <- sapply(cor_ML_sen, function(x) x$p.value)

# combine the results into a data frame
cor__ML_sendf <- data.frame(gene_id = names(cor_ML_sen), correlation = cor_coeffs_ML_sen, adj.p.value = p_values_ML_sen)

# remove the ".rho" characters from the row names
rownames(cor__ML_sendf) <- sub("\\.rho$", "", rownames(cor__ML_sendf))

hist(cor__ML_sendf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Senescence gene counts in ML")

# order the rows based on the correlation values
cor__ML_sendf <- cor__ML_sendf[order(cor__ML_sendf$correlation, decreasing = TRUE),]

fn_out_2 <- file.path(dir_outputs, "Senescence_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_sendf,fn_out_2, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Calculat Spearman rank corr
cor_CA3_sen <- apply(logcounts(spe_CA3_sen), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_sen <- sapply(cor_CA3_sen, function(x) x$estimate)
p_values_CA3_sen <- sapply(cor_CA3_sen, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_sendf <- data.frame(gene_id = names(cor_CA3_sen), correlation = cor_coeffs_CA3_sen, adj.p.value = p_values_CA3_sen)

# remove the ".rho" characters from the row names
rownames(cor__CA3_sendf) <- sub("\\.rho$", "", rownames(cor__CA3_sendf))

hist(cor__CA3_sendf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Senescence gene counts in CA3/4")

# order the rows based on the correlation values
cor__CA3_sendf <- cor__CA3_sendf[order(cor__CA3_sendf$correlation, decreasing = TRUE),]

fn_out_3 <- file.path(dir_outputs, "Senescence_CA3&4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_sendf,fn_out_3, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Calculat Spearman rank corr
cor_SGZ_sen <- apply(logcounts(spe_SGZ_sen), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_sen <- sapply(cor_SGZ_sen, function(x) x$estimate)
p_values_SGZ_sen <- sapply(cor_SGZ_sen, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_sendf <- data.frame(gene_id = names(cor_SGZ_sen), correlation = cor_coeffs_SGZ_sen, adj.p.value = p_values_SGZ_sen)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_sendf) <- sub("\\.rho$", "", rownames(cor__SGZ_sendf))

hist(cor__SGZ_sendf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Senescence gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_sendf <- cor__SGZ_sendf[order(cor__SGZ_sendf$correlation, decreasing = TRUE),]

fn_out_4 <- file.path(dir_outputs, "Senescence_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_sendf,fn_out_4, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Calculat Spearman rank corr
cor_GCL_sen <- apply(logcounts(spe_GCL_sen), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_sen <- sapply(cor_GCL_sen, function(x) x$estimate)
p_values_GCL_sen <- sapply(cor_GCL_sen, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_sendf <- data.frame(gene_id = names(cor_GCL_sen), correlation = cor_coeffs_GCL_sen, adj.p.value = p_values_GCL_sen)

# remove the ".rho" characters from the row names
rownames(cor__GCL_sendf) <- sub("\\.rho$", "", rownames(cor__GCL_sendf))

hist(cor__GCL_sendf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Senescence gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_sendf <- cor__GCL_sendf[order(cor__GCL_sendf$correlation, decreasing = TRUE),]

fn_out_5 <- file.path(dir_outputs, "Senescence_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_sendf,fn_out_5, row.names = TRUE)
