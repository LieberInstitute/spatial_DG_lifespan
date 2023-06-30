##################################################################
# spatial_DG_lifespan project
# Spearman rank corr for gene expression vs. age for SLM WM, & DG
# Anthony Ramnauth, May 22 2023
##################################################################

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

spe_WM <- spe[, which(spe$bayesSpace_harmony_10 == "10")]
dim(spe_WM)

spe_ML <- spe[, which(spe$bayesSpace_harmony_10 == "2")]
dim(spe_ML)

spe_CA3_4 <- spe[, which(spe$bayesSpace_harmony_10 == "4")]
dim(spe_CA3_4)

spe_SGZ <- spe[, which(spe$bayesSpace_harmony_10 == "6")]
dim(spe_SGZ)

spe_GCL <- spe[, which(spe$bayesSpace_harmony_10 == "7")]
dim(spe_GCL)

# Get list of gene-set from mouse data (Hansruedi Mathys et al., 2017) for microglia states

Mathys_2017 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Mathys_2017.csv"))

Su_microg_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Su_microg_2022.csv"))

Su_microg_2022 <- Su_microg_2022[Su_microg_2022$Cluster.ID == "MG1",]

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
clust_3_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.3, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_7_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.7, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_6_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.6, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Su_list <- bitr(Su_microg_2022$Gene, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

#SLM
ages_SLM <- spe_SLM$age

#WM
ages_WM <- spe_WM$age

#ML
ages_ML <- spe_ML$age

#CA3&4
ages_CA3 <- spe_CA3_4$age

#SGZ
ages_SGZ <- spe_SGZ$age

#GCL
ages_GCL <- spe_GCL$age

# Subset for Microglia clusters
spe_SLM_cl3 <- spe_SLM[which((rownames(spe_SLM)) %in% clust_3_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_SLM_cl3) <- rowData(spe_SLM_cl3)$gene_name

spe_SLM_cl6 <- spe_SLM[which((rownames(spe_SLM)) %in% clust_6_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_SLM_cl6) <- rowData(spe_SLM_cl6)$gene_name

spe_SLM_cl7 <- spe_SLM[which((rownames(spe_SLM)) %in% clust_7_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_SLM_cl7) <- rowData(spe_SLM_cl7)$gene_name

spe_WM_cl3 <- spe_WM[which((rownames(spe_WM)) %in% clust_3_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_WM_cl3) <- rowData(spe_WM_cl3)$gene_name

spe_WM_cl6 <- spe_WM[which((rownames(spe_WM)) %in% clust_6_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_WM_cl6) <- rowData(spe_WM_cl6)$gene_name

spe_WM_cl7 <- spe_WM[which((rownames(spe_WM)) %in% clust_7_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_WM_cl7) <- rowData(spe_WM_cl7)$gene_name

spe_ML_cl3 <- spe_ML[which((rownames(spe_ML)) %in% clust_3_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl3) <- rowData(spe_ML_cl3)$gene_name

spe_ML_cl6 <- spe_ML[which((rownames(spe_ML)) %in% clust_6_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl6) <- rowData(spe_ML_cl6)$gene_name

spe_ML_cl7 <- spe_ML[which((rownames(spe_ML)) %in% clust_7_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl7) <- rowData(spe_ML_cl7)$gene_name

spe_ML_Su <- spe_ML[which((rownames(spe_ML)) %in% Su_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_Su) <- rowData(spe_ML_Su)$gene_name

spe_CA3_cl3 <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% clust_3_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_CA3_cl3) <- rowData(spe_CA3_cl3)$gene_name

spe_CA3_cl6 <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% clust_6_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_CA3_cl6) <- rowData(spe_CA3_cl6)$gene_name

spe_CA3_cl7 <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% clust_7_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_CA3_cl7) <- rowData(spe_CA3_cl7)$gene_name

spe_SGZ_cl3 <- spe_SGZ[which((rownames(spe_SGZ)) %in% clust_3_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_SGZ_cl3) <- rowData(spe_SGZ_cl3)$gene_name

spe_SGZ_cl6 <- spe_SGZ[which((rownames(spe_SGZ)) %in% clust_6_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_SGZ_cl6) <- rowData(spe_SGZ_cl6)$gene_name

spe_SGZ_cl7 <- spe_SGZ[which((rownames(spe_SGZ)) %in% clust_7_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_SGZ_cl7) <- rowData(spe_SGZ_cl7)$gene_name

spe_GCL_cl3 <- spe_GCL[which((rownames(spe_GCL)) %in% clust_3_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_GCL_cl3) <- rowData(spe_GCL_cl3)$gene_name

spe_GCL_cl6 <- spe_GCL[which((rownames(spe_GCL)) %in% clust_6_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_GCL_cl6) <- rowData(spe_GCL_cl6)$gene_name

spe_GCL_cl7 <- spe_GCL[which((rownames(spe_GCL)) %in% clust_7_list$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_GCL_cl7) <- rowData(spe_GCL_cl7)$gene_name

#Cluster 3 microglia
#Calculat Spearman rank corr
cor_SLM_cl3 <- apply(logcounts(spe_SLM_cl3), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_cl3 <- sapply(cor_SLM_cl3, function(x) x$estimate)
p_values_SLM_cl3 <- sapply(cor_SLM_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_cl3df <- data.frame(gene_id = names(cor_SLM_cl3), correlation = cor_coeffs_SLM_cl3, adj.p.value = p_values_SLM_cl3)

# remove the ".rho" characters from the row names
rownames(cor__SLM_cl3df) <- sub("\\.rho$", "", rownames(cor__SLM_cl3df))

hist(cor__SLM_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_1 Microglia gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_cl3df <- cor__SLM_cl3df[order(cor__SLM_cl3df$correlation, decreasing = TRUE),]

fn_out_1 <- file.path(dir_outputs, "Early_Activated_1_Microglia_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_cl3df,fn_out_1, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 7 microglia
#Calculat Spearman rank corr
cor_SLM_cl7 <- apply(logcounts(spe_SLM_cl7), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_cl7 <- sapply(cor_SLM_cl7, function(x) x$estimate)
p_values_SLM_cl7 <- sapply(cor_SLM_cl7, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_cl7df <- data.frame(gene_id = names(cor_SLM_cl7), correlation = cor_coeffs_SLM_cl7, adj.p.value = p_values_SLM_cl7)

# remove the ".rho" characters from the row names
rownames(cor__SLM_cl7df) <- sub("\\.rho$", "", rownames(cor__SLM_cl7df))

hist(cor__SLM_cl7df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_2 Microglia gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_cl7df <- cor__SLM_cl7df[order(cor__SLM_cl7df$correlation, decreasing = TRUE),]

fn_out_2 <- file.path(dir_outputs, "Early_Activated_2_Microglia_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_cl7df,fn_out_2, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 6 microglia
#Calculat Spearman rank corr
cor_SLM_cl6 <- apply(logcounts(spe_SLM_cl6), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_cl6 <- sapply(cor_SLM_cl6, function(x) x$estimate)
p_values_SLM_cl6 <- sapply(cor_SLM_cl6, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_cl6df <- data.frame(gene_id = names(cor_SLM_cl6), correlation = cor_coeffs_SLM_cl6, adj.p.value = p_values_SLM_cl6)

# remove the ".rho" characters from the row names
rownames(cor__SLM_cl6df) <- sub("\\.rho$", "", rownames(cor__SLM_cl6df))

hist(cor__SLM_cl6df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Late Activated Microglia gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_cl6df <- cor__SLM_cl6df[order(cor__SLM_cl6df$correlation, decreasing = TRUE),]

fn_out_3 <- file.path(dir_outputs, "Late_Activated_Microglia_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_cl6df,fn_out_3, row.names = TRUE)
#----------------------------------------------------------------------------------------------------------------------------------------
#Cluster 3 microglia
#Calculat Spearman rank corr
cor_WM_cl3 <- apply(logcounts(spe_WM_cl3), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_cl3 <- sapply(cor_WM_cl3, function(x) x$estimate)
p_values_WM_cl3 <- sapply(cor_WM_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__WM_cl3df <- data.frame(gene_id = names(cor_WM_cl3), correlation = cor_coeffs_WM_cl3, adj.p.value = p_values_WM_cl3)

# remove the ".rho" characters from the row names
rownames(cor__WM_cl3df) <- sub("\\.rho$", "", rownames(cor__WM_cl3df))

hist(cor__WM_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_1 Microglia gene counts in WM")

# order the rows based on the correlation values
cor__WM_cl3df <- cor__WM_cl3df[order(cor__WM_cl3df$correlation, decreasing = TRUE),]

fn_out_4 <- file.path(dir_outputs, "Early_Activated_1_Microglia_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_cl3df,fn_out_4, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 7 microglia
#Calculat Spearman rank corr
cor_WM_cl7 <- apply(logcounts(spe_WM_cl7), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_cl7 <- sapply(cor_WM_cl7, function(x) x$estimate)
p_values_WM_cl7 <- sapply(cor_WM_cl7, function(x) x$p.value)

# combine the results into a data frame
cor__WM_cl7df <- data.frame(gene_id = names(cor_WM_cl7), correlation = cor_coeffs_WM_cl7, adj.p.value = p_values_WM_cl7)

# remove the ".rho" characters from the row names
rownames(cor__WM_cl7df) <- sub("\\.rho$", "", rownames(cor__WM_cl7df))

hist(cor__WM_cl7df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_2 Microglia gene counts in WM")

# order the rows based on the correlation values
cor__WM_cl7df <- cor__WM_cl7df[order(cor__WM_cl7df$correlation, decreasing = TRUE),]

fn_out_5 <- file.path(dir_outputs, "Early_Activated_2_Microglia_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_cl7df,fn_out_5, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 6 microglia
#Calculat Spearman rank corr
cor_WM_cl6 <- apply(logcounts(spe_WM_cl6), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_cl6 <- sapply(cor_WM_cl6, function(x) x$estimate)
p_values_WM_cl6 <- sapply(cor_WM_cl6, function(x) x$p.value)

# combine the results into a data frame
cor__WM_cl6df <- data.frame(gene_id = names(cor_WM_cl6), correlation = cor_coeffs_WM_cl6, adj.p.value = p_values_WM_cl6)

# remove the ".rho" characters from the row names
rownames(cor__WM_cl6df) <- sub("\\.rho$", "", rownames(cor__WM_cl6df))

hist(cor__WM_cl6df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Late Activated Microglia gene counts in WM")

# order the rows based on the correlation values
cor__WM_cl6df <- cor__WM_cl6df[order(cor__WM_cl6df$correlation, decreasing = TRUE),]

fn_out_6 <- file.path(dir_outputs, "Late_Activated_Microglia_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_cl6df,fn_out_6, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 3 microglia
#Calculat Spearman rank corr
cor_ML_cl3 <- apply(logcounts(spe_ML_cl3), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl3 <- sapply(cor_ML_cl3, function(x) x$estimate)
p_values_ML_cl3 <- sapply(cor_ML_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl3df <- data.frame(gene_id = names(cor_ML_cl3), correlation = cor_coeffs_ML_cl3, adj.p.value = p_values_ML_cl3)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl3df) <- sub("\\.rho$", "", rownames(cor__ML_cl3df))

hist(cor__ML_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_1 Microglia gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl3df <- cor__ML_cl3df[order(cor__ML_cl3df$correlation, decreasing = TRUE),]

fn_out_7 <- file.path(dir_outputs, "Early_Activated_1_Microglia_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl3df,fn_out_7, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 7 microglia
#Calculat Spearman rank corr
cor_ML_cl7 <- apply(logcounts(spe_ML_cl7), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl7 <- sapply(cor_ML_cl7, function(x) x$estimate)
p_values_ML_cl7 <- sapply(cor_ML_cl7, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl7df <- data.frame(gene_id = names(cor_ML_cl7), correlation = cor_coeffs_ML_cl7, adj.p.value = p_values_ML_cl7)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl7df) <- sub("\\.rho$", "", rownames(cor__ML_cl7df))

hist(cor__ML_cl7df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_2 Microglia gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl7df <- cor__ML_cl7df[order(cor__ML_cl7df$correlation, decreasing = TRUE),]

fn_out_8 <- file.path(dir_outputs, "Early_Activated_2_Microglia_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl7df,fn_out_8, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 6 microglia
#Calculat Spearman rank corr
cor_ML_cl6 <- apply(logcounts(spe_ML_cl6), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl6 <- sapply(cor_ML_cl6, function(x) x$estimate)
p_values_ML_cl6 <- sapply(cor_ML_cl6, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl6df <- data.frame(gene_id = names(cor_ML_cl6), correlation = cor_coeffs_ML_cl6, adj.p.value = p_values_ML_cl6)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl6df) <- sub("\\.rho$", "", rownames(cor__ML_cl6df))

hist(cor__ML_cl6df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Late Activated Microglia gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl6df <- cor__ML_cl6df[order(cor__ML_cl6df$correlation, decreasing = TRUE),]

fn_out_9 <- file.path(dir_outputs, "Late_Activated_Microglia_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl6df,fn_out_9, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Su_2022 microglia
#Calculate Spearman rank corr
cor_ML_Su <- apply(logcounts(spe_ML_Su), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_Su <- sapply(cor_ML_Su, function(x) x$estimate)
p_values_ML_Su <- sapply(cor_ML_Su, function(x) x$p.value)

# combine the results into a data frame
cor__ML_Sudf <- data.frame(gene_id = names(cor_ML_Su), correlation = cor_coeffs_ML_Su, adj.p.value = p_values_ML_Su)

# remove the ".rho" characters from the row names
rownames(cor__ML_Sudf) <- sub("\\.rho$", "", rownames(cor__ML_Sudf))

hist(cor__ML_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Su et al., 2022 MG1 gene counts in ML")

# order the rows based on the correlation values
cor__ML_Sudf <- cor__ML_Sudf[order(cor__ML_Sudf$correlation, decreasing = TRUE),]

fn_out_Su <- file.path(dir_outputs, "Su_2022_Microglia_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_Sudf,fn_out_Su, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------

#Cluster 3 microglia
#Calculat Spearman rank corr
cor_CA3_cl3 <- apply(logcounts(spe_CA3_cl3), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_cl3 <- sapply(cor_CA3_cl3, function(x) x$estimate)
p_values_CA3_cl3 <- sapply(cor_CA3_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_cl3df <- data.frame(gene_id = names(cor_CA3_cl3), correlation = cor_coeffs_CA3_cl3, adj.p.value = p_values_CA3_cl3)

# remove the ".rho" characters from the row names
rownames(cor__CA3_cl3df) <- sub("\\.rho$", "", rownames(cor__CA3_cl3df))

hist(cor__CA3_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_1 Microglia gene counts in CA3/4")

# order the rows based on the correlation values
cor__CA3_cl3df <- cor__CA3_cl3df[order(cor__CA3_cl3df$correlation, decreasing = TRUE),]

fn_out_10 <- file.path(dir_outputs, "Early_Activated_1_Microglia_CA3_4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_cl3df,fn_out_10, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 7 microglia
#Calculat Spearman rank corr
cor_CA3_cl7 <- apply(logcounts(spe_CA3_cl7), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_cl7 <- sapply(cor_CA3_cl7, function(x) x$estimate)
p_values_CA3_cl7 <- sapply(cor_CA3_cl7, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_cl7df <- data.frame(gene_id = names(cor_CA3_cl7), correlation = cor_coeffs_CA3_cl7, adj.p.value = p_values_CA3_cl7)

# remove the ".rho" characters from the row names
rownames(cor__CA3_cl7df) <- sub("\\.rho$", "", rownames(cor__CA3_cl7df))

hist(cor__CA3_cl7df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_2 Microglia gene counts in CA3/4")

# order the rows based on the correlation values
cor__CA3_cl7df <- cor__CA3_cl7df[order(cor__CA3_cl7df$correlation, decreasing = TRUE),]

fn_out_11 <- file.path(dir_outputs, "Early_Activated_2_Microglia_CA3_4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_cl7df,fn_out_11, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 6 microglia
#Calculat Spearman rank corr
cor_CA3_cl6 <- apply(logcounts(spe_CA3_cl6), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_cl6 <- sapply(cor_CA3_cl6, function(x) x$estimate)
p_values_CA3_cl6 <- sapply(cor_CA3_cl6, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_cl6df <- data.frame(gene_id = names(cor_CA3_cl6), correlation = cor_coeffs_CA3_cl6, adj.p.value = p_values_CA3_cl6)

# remove the ".rho" characters from the row names
rownames(cor__CA3_cl6df) <- sub("\\.rho$", "", rownames(cor__CA3_cl6df))

hist(cor__CA3_cl6df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Late Activated Microglia gene counts in CA3/4")

# order the rows based on the correlation values
cor__CA3_cl6df <- cor__CA3_cl6df[order(cor__CA3_cl6df$correlation, decreasing = TRUE),]

fn_out_12 <- file.path(dir_outputs, "Late_Activated_Microglia_CA3_4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_cl6df,fn_out_12, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 3 microglia
#Calculat Spearman rank corr
cor_SGZ_cl3 <- apply(logcounts(spe_SGZ_cl3), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_cl3 <- sapply(cor_SGZ_cl3, function(x) x$estimate)
p_values_SGZ_cl3 <- sapply(cor_SGZ_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_cl3df <- data.frame(gene_id = names(cor_SGZ_cl3), correlation = cor_coeffs_SGZ_cl3, adj.p.value = p_values_SGZ_cl3)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_cl3df) <- sub("\\.rho$", "", rownames(cor__SGZ_cl3df))

hist(cor__SGZ_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_1 Microglia gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_cl3df <- cor__SGZ_cl3df[order(cor__SGZ_cl3df$correlation, decreasing = TRUE),]

fn_out_13 <- file.path(dir_outputs, "Early_Activated_1_Microglia_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_cl3df,fn_out_13, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 7 microglia
#Calculat Spearman rank corr
cor_SGZ_cl7 <- apply(logcounts(spe_SGZ_cl7), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_cl7 <- sapply(cor_SGZ_cl7, function(x) x$estimate)
p_values_SGZ_cl7 <- sapply(cor_SGZ_cl7, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_cl7df <- data.frame(gene_id = names(cor_SGZ_cl7), correlation = cor_coeffs_SGZ_cl7, adj.p.value = p_values_SGZ_cl7)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_cl7df) <- sub("\\.rho$", "", rownames(cor__SGZ_cl7df))

hist(cor__SGZ_cl7df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_2 Microglia gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_cl7df <- cor__SGZ_cl7df[order(cor__SGZ_cl7df$correlation, decreasing = TRUE),]

fn_out_14 <- file.path(dir_outputs, "Early_Activated_2_Microglia_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_cl7df,fn_out_14, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 6 microglia
#Calculat Spearman rank corr
cor_SGZ_cl6 <- apply(logcounts(spe_SGZ_cl6), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_cl6 <- sapply(cor_SGZ_cl6, function(x) x$estimate)
p_values_SGZ_cl6 <- sapply(cor_SGZ_cl6, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_cl6df <- data.frame(gene_id = names(cor_SGZ_cl6), correlation = cor_coeffs_SGZ_cl6, adj.p.value = p_values_SGZ_cl6)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_cl6df) <- sub("\\.rho$", "", rownames(cor__SGZ_cl6df))

hist(cor__SGZ_cl6df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Late Activated Microglia gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_cl6df <- cor__SGZ_cl6df[order(cor__SGZ_cl6df$correlation, decreasing = TRUE),]

fn_out_15 <- file.path(dir_outputs, "Late_Activated_Microglia_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_cl6df,fn_out_15, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 3 microglia
#Calculat Spearman rank corr
cor_GCL_cl3 <- apply(logcounts(spe_GCL_cl3), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_cl3 <- sapply(cor_GCL_cl3, function(x) x$estimate)
p_values_GCL_cl3 <- sapply(cor_GCL_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_cl3df <- data.frame(gene_id = names(cor_GCL_cl3), correlation = cor_coeffs_GCL_cl3, adj.p.value = p_values_GCL_cl3)

# remove the ".rho" characters from the row names
rownames(cor__GCL_cl3df) <- sub("\\.rho$", "", rownames(cor__GCL_cl3df))

hist(cor__GCL_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_1 Microglia gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_cl3df <- cor__GCL_cl3df[order(cor__GCL_cl3df$correlation, decreasing = TRUE),]

fn_out_16 <- file.path(dir_outputs, "Early_Activated_1_Microglia_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_cl3df,fn_out_16, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 7 microglia
#Calculat Spearman rank corr
cor_GCL_cl7 <- apply(logcounts(spe_GCL_cl7), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_cl7 <- sapply(cor_GCL_cl7, function(x) x$estimate)
p_values_GCL_cl7 <- sapply(cor_GCL_cl7, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_cl7df <- data.frame(gene_id = names(cor_GCL_cl7), correlation = cor_coeffs_GCL_cl7, adj.p.value = p_values_GCL_cl7)

# remove the ".rho" characters from the row names
rownames(cor__GCL_cl7df) <- sub("\\.rho$", "", rownames(cor__GCL_cl7df))

hist(cor__GCL_cl7df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Early Activated_2 Microglia gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_cl7df <- cor__GCL_cl7df[order(cor__GCL_cl7df$correlation, decreasing = TRUE),]

fn_out_17 <- file.path(dir_outputs, "Early_Activated_2_Microglia_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_cl7df,fn_out_17, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 6 microglia
#Calculat Spearman rank corr
cor_GCL_cl6 <- apply(logcounts(spe_GCL_cl6), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_cl6 <- sapply(cor_GCL_cl6, function(x) x$estimate)
p_values_GCL_cl6 <- sapply(cor_GCL_cl6, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_cl6df <- data.frame(gene_id = names(cor_GCL_cl6), correlation = cor_coeffs_GCL_cl6, adj.p.value = p_values_GCL_cl6)

# remove the ".rho" characters from the row names
rownames(cor__GCL_cl6df) <- sub("\\.rho$", "", rownames(cor__GCL_cl6df))

hist(cor__GCL_cl6df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Late Activated Microglia gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_cl6df <- cor__GCL_cl6df[order(cor__GCL_cl6df$correlation, decreasing = TRUE),]

fn_out_18 <- file.path(dir_outputs, "Late_Activated_Microglia_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_cl6df,fn_out_18, row.names = TRUE)
########################################################################################################################################
