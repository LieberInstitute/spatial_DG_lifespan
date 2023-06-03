##################################################################
# spatial_DG_lifespan project
# Spearman rank corr for dendritic gene expression vs. age in ML
# Anthony Ramnauth, May 30 2023
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

spe_ML <- spe[, which(spe$bayesSpace_harmony_10 == "2")]
dim(spe_ML)

# Get list of gene-set from mouse data (Robert R. Stickels et al., 2020) for dendritically enriched gene sets
Stickels_2020 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Stickels_2020.csv"))

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Isolated gene set from looking at radial distances from GCL boarder to ML
combined_gene_set <- c(
    "CAMK2A", "TUBB2B", "EEF1A1", "TULP4", "LHX2", "NCS1", "NSMF", "GLUD1", "PSD",
    "RDX", "DDN", "ARHGAP5", "HOMER2", "RPL13", "RPL26", "PPP1R9B", "RPS15", "EEF2",
    "SLC25A23", "SIRT2", "DACT3", "RPL13A", "CENPB", "CABP7", "RPL22", "MAP4K4",
    "MAP2", "CCNI", "RPS3A", "RPS14", "RPS24", "DENND5A", "MAPK8IP1", "FTH1", "EEF1G",
    "COX8A", "PITPNM1", "RPS3", "AGAP2", "RPLP0", "COX6A1", "CEP170B", "RPLP1", "FBXL16",
    "RPS2", "COX4I1", "RPL23", "ABHD17A", "RPS11", "RPS9", "RPS5", "ITSN1", "GPM6B",
    "IDS"
)

# Translate from one species to the other using the orthology
Dcluster_1 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.1,]
Dcluster_2 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.2,]
Dcluster_3 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.3,]
Dcluster_4 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster4,]

#Get the Ensembl IDs
Dcluster_1 <- bitr(Dcluster_1$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_2 <- bitr(Dcluster_2$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_3 <- bitr(Dcluster_3$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_4 <- bitr(Dcluster_4$Column1, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)
Dcluster_c <- bitr(combined_gene_set, fromType="SYMBOL", toType = "ENSEMBL", OrgDb=org.Hs.eg.db)

#ML age
ages_ML <- spe_ML$age

# Subset for dendritic clusters
spe_ML_cl1 <- spe_ML[which((rownames(spe_ML)) %in% Dcluster_1$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl1) <- rowData(spe_ML_cl1)$gene_name

spe_ML_cl2 <- spe_ML[which((rownames(spe_ML)) %in% Dcluster_2$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl2) <- rowData(spe_ML_cl2)$gene_name

spe_ML_cl3 <- spe_ML[which((rownames(spe_ML)) %in% Dcluster_3$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl3) <- rowData(spe_ML_cl3)$gene_name

spe_ML_cl4 <- spe_ML[which((rownames(spe_ML)) %in% Dcluster_4$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_cl4) <- rowData(spe_ML_cl4)$gene_name

# Subset for dendritic clusters
spe_ML_clc <- spe_ML[which((rownames(spe_ML)) %in% Dcluster_c$ENSEMBL)]
# Set gene names as row names for easier plotting
rownames(spe_ML_clc) <- rowData(spe_ML_clc)$gene_name

########################################################################################################################################

#Cluster 1 dendritic gene set
#Calculat Spearman rank corr
cor_ML_cl1 <- apply(logcounts(spe_ML_cl1), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl1 <- sapply(cor_ML_cl1, function(x) x$estimate)
p_values_ML_cl1 <- sapply(cor_ML_cl1, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl1df <- data.frame(gene_id = names(cor_ML_cl1), correlation = cor_coeffs_ML_cl1, adj.p.value = p_values_ML_cl1)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl1df) <- sub("\\.rho$", "", rownames(cor__ML_cl1df))

hist(cor__ML_cl1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Dendritic enriched cluster 1 gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl1df <- cor__ML_cl1df[order(cor__ML_cl1df$correlation, decreasing = TRUE),]

fn_out_1 <- file.path(dir_outputs, "Dendritic_enriched1_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl1df,fn_out_1, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 2 dendritic gene set
#Calculat Spearman rank corr
cor_ML_cl2 <- apply(logcounts(spe_ML_cl2), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl2 <- sapply(cor_ML_cl2, function(x) x$estimate)
p_values_ML_cl2 <- sapply(cor_ML_cl2, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl2df <- data.frame(gene_id = names(cor_ML_cl2), correlation = cor_coeffs_ML_cl2, adj.p.value = p_values_ML_cl2)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl2df) <- sub("\\.rho$", "", rownames(cor__ML_cl2df))

hist(cor__ML_cl2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Dendritic enriched cluster 2 gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl2df <- cor__ML_cl2df[order(cor__ML_cl2df$correlation, decreasing = TRUE),]

fn_out_2 <- file.path(dir_outputs, "Dendritic_enriched2_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl2df,fn_out_2, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 3 dendritic gene set
#Calculat Spearman rank corr
cor_ML_cl3 <- apply(logcounts(spe_ML_cl3), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl3 <- sapply(cor_ML_cl3, function(x) x$estimate)
p_values_ML_cl3 <- sapply(cor_ML_cl3, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl3df <- data.frame(gene_id = names(cor_ML_cl3), correlation = cor_coeffs_ML_cl3, adj.p.value = p_values_ML_cl3)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl3df) <- sub("\\.rho$", "", rownames(cor__ML_cl3df))

hist(cor__ML_cl3df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Dendritic enriched cluster 3 gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl3df <- cor__ML_cl3df[order(cor__ML_cl3df$correlation, decreasing = TRUE),]

fn_out_3 <- file.path(dir_outputs, "Dendritic_enriched3_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl3df,fn_out_3, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Cluster 4 dendritic gene set
#Calculat Spearman rank corr
cor_ML_cl4 <- apply(logcounts(spe_ML_cl4), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_cl4 <- sapply(cor_ML_cl4, function(x) x$estimate)
p_values_ML_cl4 <- sapply(cor_ML_cl4, function(x) x$p.value)

# combine the results into a data frame
cor__ML_cl4df <- data.frame(gene_id = names(cor_ML_cl4), correlation = cor_coeffs_ML_cl4, adj.p.value = p_values_ML_cl4)

# remove the ".rho" characters from the row names
rownames(cor__ML_cl4df) <- sub("\\.rho$", "", rownames(cor__ML_cl4df))

hist(cor__ML_cl4df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Dendritic enriched cluster 4 gene counts in ML")

# order the rows based on the correlation values
cor__ML_cl4df <- cor__ML_cl4df[order(cor__ML_cl4df$correlation, decreasing = TRUE),]

fn_out_4 <- file.path(dir_outputs, "Dendritic_enriched4_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_cl4df,fn_out_4, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------

# Combined radial dist. from GCL verified dendritic gene set
#Calculat Spearman rank corr
cor_ML_clc <- apply(logcounts(spe_ML_clc), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_clc <- sapply(cor_ML_clc, function(x) x$estimate)
p_values_ML_clc <- sapply(cor_ML_clc, function(x) x$p.value)

# combine the results into a data frame
cor__ML_clcdf <- data.frame(gene_id = names(cor_ML_clc), correlation = cor_coeffs_ML_clc, adj.p.value = p_values_ML_clc)

# remove the ".rho" characters from the row names
rownames(cor__ML_clcdf) <- sub("\\.rho$", "", rownames(cor__ML_clcdf))

hist(cor__ML_clcdf$correlation, xlab = "age rho correlation", ylab = "Gene count",
    main = "Dendritic enriched radial dist. away from GCL gene counts in ML")

# order the rows based on the correlation values
cor__ML_clcdf <- cor__ML_clcdf[order(cor__ML_clcdf$correlation, decreasing = TRUE),]

fn_out_c <- file.path(dir_outputs, "Dendritic_enriched_radial_dist_GCL_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_clcdf,fn_out_c, row.names = TRUE)

