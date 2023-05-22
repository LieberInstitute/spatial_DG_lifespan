##################################################################
# spatial_DG_lifespan project
# Spearman rank corr for gene expression vs. age for each SLM & WM
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

# Get list of gene-set from mouse data (Hansruedi Mathys et al., 2017) for microglia states

Mathys_2017 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Mathys_2017.csv"))

# Convert to Ensembl IDs for gene_set_enrichment funciton to work
clust_3_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.3, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_7_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.7, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
clust_6_list <- bitr(Mathys_2017$Genes.up.regulated.in.Cluster.6, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

#ML
ages_SLM <- spe_SLM$age

#CA3&4
ages_WM <- spe_WM$age

# Subset for Microglia clusters
spe_SLM_cl3 <- spe_SLM[which((rownames(spe_SLM)) %in% clust_3_list$ENSEMBL)]
spe_SLM_cl6 <- spe_SLM[which((rownames(spe_SLM)) %in% clust_6_list$ENSEMBL)]
spe_SLM_cl7 <- spe_SLM[which((rownames(spe_SLM)) %in% clust_7_list$ENSEMBL)]
spe_WM_cl3 <- spe_WM[which((rownames(spe_WM)) %in% clust_3_list$ENSEMBL)]
spe_WM_cl6 <- spe_WM[which((rownames(spe_WM)) %in% clust_6_list$ENSEMBL)]
spe_WM_cl7 <- spe_WM[which((rownames(spe_WM)) %in% clust_7_list$ENSEMBL)]

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

