###########################################################################
# spatial_DG_lifespan project
# Spearman rank corr for gene expression vs. age for SLM WM, & DG for astro
# Anthony Ramnauth, June 22 2023
###########################################################################

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

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

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

# Reactive Astrocyte gene sets

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Manually input Pan, A1, & A2 gene sets from mouse data (Laura Clarke et al., 2018)
Clarke_2018 <- list(
    PAN = c("Lcn2", "Steap4", "S1pr3", "Timp1", "Hsbp1", "Cxcl10", "Cd44", "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap"),
    A1 = c("C3", "H2-T23", "Serping1", "H2-D1", "Ggta1", "Ligp1", "Gpp2", "Fbln5", "Fkbp5", "Psmb8", "Srgn", "Amigo2"),
    A2 = c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14", "Stat3")
)

# Translate from one species to the other using the orthology
PAN <- orthology[orthology$Column3 %in% Clarke_2018$PAN,]
A1 <- orthology[orthology$Column3 %in% Clarke_2018$A1,]
A2 <- orthology[orthology$Column3 %in% Clarke_2018$A2,]

# Get list of gene-set from human data (Yijing Su et al., 2022) for Reactive astrocytes
Su_astro_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Su_astro_2022.csv"))

Su_astro1_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST1",]
Su_astro6_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST6",]
Su_astro_2022 <- rbind(Su_astro1_2022, Su_astro6_2022)

# Subset for Astrocyte clusters
spe_SLM_PAN <- spe_SLM[which((rownames(spe_SLM)) %in% PAN$Column1)]
spe_SLM_A1 <- spe_SLM[which((rownames(spe_SLM)) %in% A1$Column1)]
spe_SLM_A2 <- spe_SLM[which((rownames(spe_SLM)) %in% A2$Column1)]
spe_SLM_Su <- spe_SLM[which((rownames(spe_SLM)) %in% Su_astro_2022$Gene)]
spe_WM_PAN <- spe_WM[which((rownames(spe_WM)) %in% PAN$Column1)]
spe_WM_A1 <- spe_WM[which((rownames(spe_WM)) %in% A1$Column1)]
spe_WM_A2 <- spe_WM[which((rownames(spe_WM)) %in% A2$Column1)]
spe_WM_Su <- spe_WM[which((rownames(spe_WM)) %in% Su_astro_2022$Gene)]
spe_ML_PAN <- spe_ML[which((rownames(spe_ML)) %in% PAN$Column1)]
spe_ML_A1 <- spe_ML[which((rownames(spe_ML)) %in% A1$Column1)]
spe_ML_A2 <- spe_ML[which((rownames(spe_ML)) %in% A2$Column1)]
spe_ML_Su <- spe_ML[which((rownames(spe_ML)) %in% Su_astro_2022$Gene)]
spe_CA3_PAN <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% PAN$Column1)]
spe_CA3_A1 <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% A1$Column1)]
spe_CA3_A2 <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% A2$Column1)]
spe_CA3_Su <- spe_CA3_4[which((rownames(spe_CA3_4)) %in% Su_astro_2022$Gene)]
spe_SGZ_PAN <- spe_SGZ[which((rownames(spe_SGZ)) %in% PAN$Column1)]
spe_SGZ_A1 <- spe_SGZ[which((rownames(spe_SGZ)) %in% A1$Column1)]
spe_SGZ_A2 <- spe_SGZ[which((rownames(spe_SGZ)) %in% A2$Column1)]
spe_SGZ_Su <- spe_SGZ[which((rownames(spe_SGZ)) %in% Su_astro_2022$Gene)]
spe_GCL_PAN <- spe_GCL[which((rownames(spe_GCL)) %in% PAN$Column1)]
spe_GCL_A1 <- spe_GCL[which((rownames(spe_GCL)) %in% A1$Column1)]
spe_GCL_A2 <- spe_GCL[which((rownames(spe_GCL)) %in% A2$Column1)]
spe_GCL_Su <- spe_GCL[which((rownames(spe_GCL)) %in% Su_astro_2022$Gene)]

#PAN astrocyte markers
#Calculat Spearman rank corr
cor_SLM_PAN <- apply(logcounts(spe_SLM_PAN), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_PAN <- sapply(cor_SLM_PAN, function(x) x$estimate)
p_values_SLM_PAN <- sapply(cor_SLM_PAN, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_PANdf <- data.frame(gene_id = names(cor_SLM_PAN), correlation = cor_coeffs_SLM_PAN, adj.p.value = p_values_SLM_PAN)

# remove the ".rho" characters from the row names
rownames(cor__SLM_PANdf) <- sub("\\.rho$", "", rownames(cor__SLM_PANdf))

hist(cor__SLM_PANdf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "PAN reactive Astrocyte gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_PANdf <- cor__SLM_PANdf[order(cor__SLM_PANdf$correlation, decreasing = TRUE),]

fn_out_1 <- file.path(dir_outputs, "PAN_reactive_Astrocyte_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_PANdf,fn_out_1, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A1 astrocyte markers
#Calculat Spearman rank corr
cor_SLM_A1 <- apply(logcounts(spe_SLM_A1), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_A1 <- sapply(cor_SLM_A1, function(x) x$estimate)
p_values_SLM_A1 <- sapply(cor_SLM_A1, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_A1df <- data.frame(gene_id = names(cor_SLM_A1), correlation = cor_coeffs_SLM_A1, adj.p.value = p_values_SLM_A1)

# remove the ".rho" characters from the row names
rownames(cor__SLM_A1df) <- sub("\\.rho$", "", rownames(cor__SLM_A1df))

hist(cor__SLM_A1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A1 reactive Astrocyte gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_A1df <- cor__SLM_A1df[order(cor__SLM_A1df$correlation, decreasing = TRUE),]

fn_out_2 <- file.path(dir_outputs, "A1_reactive_Astrocyte_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_A1df,fn_out_2, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A2 astrocyte markers
#Calculat Spearman rank corr
cor_SLM_A2 <- apply(logcounts(spe_SLM_A2), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_A2 <- sapply(cor_SLM_A2, function(x) x$estimate)
p_values_SLM_A2 <- sapply(cor_SLM_A2, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_A2df <- data.frame(gene_id = names(cor_SLM_A2), correlation = cor_coeffs_SLM_A2, adj.p.value = p_values_SLM_A2)

# remove the ".rho" characters from the row names
rownames(cor__SLM_A2df) <- sub("\\.rho$", "", rownames(cor__SLM_A2df))

hist(cor__SLM_A2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A2 reactive Astrocyte gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_A2df <- cor__SLM_A2df[order(cor__SLM_A2df$correlation, decreasing = TRUE),]

fn_out_3 <- file.path(dir_outputs, "A2_reactive_Astrocyte_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_A2df,fn_out_3, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Ast1&6 astrocyte markers
#Calculat Spearman rank corr
cor_SLM_Su <- apply(logcounts(spe_SLM_Su), 1, function(row) cor.test(as.numeric(row), ages_SLM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SLM_Su <- sapply(cor_SLM_Su, function(x) x$estimate)
p_values_SLM_Su <- sapply(cor_SLM_Su, function(x) x$p.value)

# combine the results into a data frame
cor__SLM_Sudf <- data.frame(gene_id = names(cor_SLM_Su), correlation = cor_coeffs_SLM_Su, adj.p.value = p_values_SLM_Su)

# remove the ".rho" characters from the row names
rownames(cor__SLM_Sudf) <- sub("\\.rho$", "", rownames(cor__SLM_Sudf))

hist(cor__SLM_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Ast1&6 reactive Astrocyte gene counts in SLM")

# order the rows based on the correlation values
cor__SLM_Sudf <- cor__SLM_Sudf[order(cor__SLM_Sudf$correlation, decreasing = TRUE),]

fn_out_4 <- file.path(dir_outputs, "Ast1&6_reactive_Astrocyte_SLM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SLM_Sudf,fn_out_4, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#PAN astrocyte markers
#Calculat Spearman rank corr
cor_WM_PAN <- apply(logcounts(spe_WM_PAN), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_PAN <- sapply(cor_WM_PAN, function(x) x$estimate)
p_values_WM_PAN <- sapply(cor_WM_PAN, function(x) x$p.value)

# combine the results into a data frame
cor__WM_PANdf <- data.frame(gene_id = names(cor_WM_PAN), correlation = cor_coeffs_WM_PAN, adj.p.value = p_values_WM_PAN)

# remove the ".rho" characters from the row names
rownames(cor__WM_PANdf) <- sub("\\.rho$", "", rownames(cor__WM_PANdf))

hist(cor__WM_PANdf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "PAN reactive Astrocyte gene counts in WM")

# order the rows based on the correlation values
cor__WM_PANdf <- cor__WM_PANdf[order(cor__WM_PANdf$correlation, decreasing = TRUE),]

fn_out_5 <- file.path(dir_outputs, "PAN_reactive_Astrocyte_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_PANdf,fn_out_5, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A1 astrocyte markers
#Calculat Spearman rank corr
cor_WM_A1 <- apply(logcounts(spe_WM_A1), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_A1 <- sapply(cor_WM_A1, function(x) x$estimate)
p_values_WM_A1 <- sapply(cor_WM_A1, function(x) x$p.value)

# combine the results into a data frame
cor__WM_A1df <- data.frame(gene_id = names(cor_WM_A1), correlation = cor_coeffs_WM_A1, adj.p.value = p_values_WM_A1)

# remove the ".rho" characters from the row names
rownames(cor__WM_A1df) <- sub("\\.rho$", "", rownames(cor__WM_A1df))

hist(cor__WM_A1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A1 reactive Astrocyte gene counts in WM")

# order the rows based on the correlation values
cor__WM_A1df <- cor__WM_A1df[order(cor__WM_A1df$correlation, decreasing = TRUE),]

fn_out_6 <- file.path(dir_outputs, "A1_reactive_Astrocyte_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_A1df,fn_out_6, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A2 astrocyte markers
#Calculat Spearman rank corr
cor_WM_A2 <- apply(logcounts(spe_WM_A2), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_A2 <- sapply(cor_WM_A2, function(x) x$estimate)
p_values_WM_A2 <- sapply(cor_WM_A2, function(x) x$p.value)

# combine the results into a data frame
cor__WM_A2df <- data.frame(gene_id = names(cor_WM_A2), correlation = cor_coeffs_WM_A2, adj.p.value = p_values_WM_A2)

# remove the ".rho" characters from the row names
rownames(cor__WM_A2df) <- sub("\\.rho$", "", rownames(cor__WM_A2df))

hist(cor__WM_A2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A2 reactive Astrocyte gene counts in WM")

# order the rows based on the correlation values
cor__WM_A2df <- cor__WM_A2df[order(cor__WM_A2df$correlation, decreasing = TRUE),]

fn_out_7 <- file.path(dir_outputs, "A2_reactive_Astrocyte_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_A2df,fn_out_7, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Ast1&6 astrocyte markers
#Calculat Spearman rank corr
cor_WM_Su <- apply(logcounts(spe_WM_Su), 1, function(row) cor.test(as.numeric(row), ages_WM, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_WM_Su <- sapply(cor_WM_Su, function(x) x$estimate)
p_values_WM_Su <- sapply(cor_WM_Su, function(x) x$p.value)

# combine the results into a data frame
cor__WM_Sudf <- data.frame(gene_id = names(cor_WM_Su), correlation = cor_coeffs_WM_Su, adj.p.value = p_values_WM_Su)

# remove the ".rho" characters from the row names
rownames(cor__WM_Sudf) <- sub("\\.rho$", "", rownames(cor__WM_Sudf))

hist(cor__WM_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Ast1&6 reactive Astrocyte gene counts in WM")

# order the rows based on the correlation values
cor__WM_Sudf <- cor__WM_Sudf[order(cor__WM_Sudf$correlation, decreasing = TRUE),]

fn_out_8 <- file.path(dir_outputs, "Ast1&6_reactive_Astrocyte_WM_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__WM_Sudf,fn_out_8, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#PAN astrocyte markers
#Calculat Spearman rank corr
cor_ML_PAN <- apply(logcounts(spe_ML_PAN), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_PAN <- sapply(cor_ML_PAN, function(x) x$estimate)
p_values_ML_PAN <- sapply(cor_ML_PAN, function(x) x$p.value)

# combine the results into a data frame
cor__ML_PANdf <- data.frame(gene_id = names(cor_SLM_PAN), correlation = cor_coeffs_ML_PAN, adj.p.value = p_values_ML_PAN)

# remove the ".rho" characters from the row names
rownames(cor__ML_PANdf) <- sub("\\.rho$", "", rownames(cor__ML_PANdf))

hist(cor__ML_PANdf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "PAN reactive Astrocyte gene counts in ML")

# order the rows based on the correlation values
cor__ML_PANdf <- cor__ML_PANdf[order(cor__ML_PANdf$correlation, decreasing = TRUE),]

fn_out_9 <- file.path(dir_outputs, "PAN_reactive_Astrocyte_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_PANdf,fn_out_9, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A1 astrocyte markers
#Calculat Spearman rank corr
cor_ML_A1 <- apply(logcounts(spe_ML_A1), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_A1 <- sapply(cor_ML_A1, function(x) x$estimate)
p_values_ML_A1 <- sapply(cor_ML_A1, function(x) x$p.value)

# combine the results into a data frame
cor__ML_A1df <- data.frame(gene_id = names(cor_ML_A1), correlation = cor_coeffs_ML_A1, adj.p.value = p_values_ML_A1)

# remove the ".rho" characters from the row names
rownames(cor__ML_A1df) <- sub("\\.rho$", "", rownames(cor__ML_A1df))

hist(cor__ML_A1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A1 reactive Astrocyte gene counts in ML")

# order the rows based on the correlation values
cor__ML_A1df <- cor__ML_A1df[order(cor__ML_A1df$correlation, decreasing = TRUE),]

fn_out_10 <- file.path(dir_outputs, "A1_reactive_Astrocyte_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_A1df,fn_out_10, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A2 astrocyte markers
#Calculat Spearman rank corr
cor_ML_A2 <- apply(logcounts(spe_ML_A2), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_A2 <- sapply(cor_ML_A2, function(x) x$estimate)
p_values_ML_A2 <- sapply(cor_ML_A2, function(x) x$p.value)

# combine the results into a data frame
cor__ML_A2df <- data.frame(gene_id = names(cor_ML_A2), correlation = cor_coeffs_ML_A2, adj.p.value = p_values_ML_A2)

# remove the ".rho" characters from the row names
rownames(cor__ML_A2df) <- sub("\\.rho$", "", rownames(cor__ML_A2df))

hist(cor__ML_A2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A2 reactive Astrocyte gene counts in ML")

# order the rows based on the correlation values
cor__ML_A2df <- cor__ML_A2df[order(cor__ML_A2df$correlation, decreasing = TRUE),]

fn_out_11 <- file.path(dir_outputs, "A2_reactive_Astrocyte_cor__ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_A2df,fn_out_11, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Ast1&6 astrocyte markers
#Calculat Spearman rank corr
cor_ML_Su <- apply(logcounts(spe_ML_Su), 1, function(row) cor.test(as.numeric(row), ages_ML, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_ML_Su <- sapply(cor_ML_Su, function(x) x$estimate)
p_values_ML_Su <- sapply(cor_ML_Su, function(x) x$p.value)

# combine the results into a data frame
cor__ML_Sudf <- data.frame(gene_id = names(cor_ML_Su), correlation = cor_coeffs_ML_Su, adj.p.value = p_values_ML_Su)

# remove the ".rho" characters from the row names
rownames(cor__ML_Sudf) <- sub("\\.rho$", "", rownames(cor__ML_Sudf))

hist(cor__ML_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Ast1&6 reactive Astrocyte gene counts in ML")

# order the rows based on the correlation values
cor__ML_Sudf <- cor__ML_Sudf[order(cor__ML_Sudf$correlation, decreasing = TRUE),]

fn_out_12 <- file.path(dir_outputs, "Ast1&6_reactive_Astrocyte_ML_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__ML_Sudf,fn_out_12, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#PAN astrocyte markers
#Calculat Spearman rank corr
cor_CA3_PAN <- apply(logcounts(spe_CA3_PAN), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_PAN <- sapply(cor_CA3_PAN, function(x) x$estimate)
p_values_CA3_PAN <- sapply(cor_CA3_PAN, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_PANdf <- data.frame(gene_id = names(cor_CA3_PAN), correlation = cor_coeffs_CA3_PAN, adj.p.value = p_values_CA3_PAN)

# remove the ".rho" characters from the row names
rownames(cor__CA3_PANdf) <- sub("\\.rho$", "", rownames(cor__CA3_PANdf))

hist(cor__CA3_PANdf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "PAN reactive Astrocyte gene counts in CA3&4")

# order the rows based on the correlation values
cor__CA3_PANdf <- cor__CA3_PANdf[order(cor__CA3_PANdf$correlation, decreasing = TRUE),]

fn_out_13 <- file.path(dir_outputs, "PAN_reactive_Astrocyte_CA3&4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_PANdf,fn_out_13, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A1 astrocyte markers
#Calculat Spearman rank corr
cor_CA3_A1 <- apply(logcounts(spe_CA3_A1), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_A1 <- sapply(cor_CA3_A1, function(x) x$estimate)
p_values_CA3_A1 <- sapply(cor_CA3_A1, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_A1df <- data.frame(gene_id = names(cor_CA3_A1), correlation = cor_coeffs_CA3_A1, adj.p.value = p_values_CA3_A1)

# remove the ".rho" characters from the row names
rownames(cor__CA3_A1df) <- sub("\\.rho$", "", rownames(cor__CA3_A1df))

hist(cor__CA3_A1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A1 reactive Astrocyte gene counts in CA3&4")

# order the rows based on the correlation values
cor__CA3_A1df <- cor__CA3_A1df[order(cor__CA3_A1df$correlation, decreasing = TRUE),]

fn_out_14 <- file.path(dir_outputs, "A1_reactive_Astrocyte_CA3&4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_A1df,fn_out_14, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A2 astrocyte markers
#Calculat Spearman rank corr
cor_CA3_A2 <- apply(logcounts(spe_CA3_A2), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_A2 <- sapply(cor_CA3_A2, function(x) x$estimate)
p_values_CA3_A2 <- sapply(cor_CA3_A2, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_A2df <- data.frame(gene_id = names(cor_ML_A2), correlation = cor_coeffs_CA3_A2, adj.p.value = p_values_CA3_A2)

# remove the ".rho" characters from the row names
rownames(cor__CA3_A2df) <- sub("\\.rho$", "", rownames(cor__CA3_A2df))

hist(cor__CA3_A2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A2 reactive Astrocyte gene counts in CA3&4")

# order the rows based on the correlation values
cor__CA3_A2df <- cor__CA3_A2df[order(cor__CA3_A2df$correlation, decreasing = TRUE),]

fn_out_15 <- file.path(dir_outputs, "A2_reactive_Astrocyte_cor__CA3&4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_A2df,fn_out_15, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Ast1&6 astrocyte markers
#Calculat Spearman rank corr
cor_CA3_Su <- apply(logcounts(spe_CA3_Su), 1, function(row) cor.test(as.numeric(row), ages_CA3, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_CA3_Su <- sapply(cor_CA3_Su, function(x) x$estimate)
p_values_CA3_Su <- sapply(cor_CA3_Su, function(x) x$p.value)

# combine the results into a data frame
cor__CA3_Sudf <- data.frame(gene_id = names(cor_CA3_Su), correlation = cor_coeffs_CA3_Su, adj.p.value = p_values_CA3_Su)

# remove the ".rho" characters from the row names
rownames(cor__CA3_Sudf) <- sub("\\.rho$", "", rownames(cor__CA3_Sudf))

hist(cor__CA3_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Ast1&6 reactive Astrocyte gene counts in CA3&4")

# order the rows based on the correlation values
cor__CA3_Sudf <- cor__CA3_Sudf[order(cor__CA3_Sudf$correlation, decreasing = TRUE),]

fn_out_16 <- file.path(dir_outputs, "Ast1&6_reactive_Astrocyte_CA3&4_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__CA3_Sudf,fn_out_16, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#PAN astrocyte markers
#Calculat Spearman rank corr
cor_SGZ_PAN <- apply(logcounts(spe_SGZ_PAN), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_PAN <- sapply(cor_SGZ_PAN, function(x) x$estimate)
p_values_SGZ_PAN <- sapply(cor_SGZ_PAN, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_PANdf <- data.frame(gene_id = names(cor_SGZ_PAN), correlation = cor_coeffs_SGZ_PAN, adj.p.value = p_values_SGZ_PAN)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_PANdf) <- sub("\\.rho$", "", rownames(cor__SGZ_PANdf))

hist(cor__SGZ_PANdf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "PAN reactive Astrocyte gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_PANdf <- cor__SGZ_PANdf[order(cor__SGZ_PANdf$correlation, decreasing = TRUE),]

fn_out_17 <- file.path(dir_outputs, "PAN_reactive_Astrocyte_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_PANdf,fn_out_17, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A1 astrocyte markers
#Calculat Spearman rank corr
cor_SGZ_A1 <- apply(logcounts(spe_SGZ_A1), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_A1 <- sapply(cor_SGZ_A1, function(x) x$estimate)
p_values_SGZ_A1 <- sapply(cor_SGZ_A1, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_A1df <- data.frame(gene_id = names(cor_SGZ_A1), correlation = cor_coeffs_SGZ_A1, adj.p.value = p_values_SGZ_A1)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_A1df) <- sub("\\.rho$", "", rownames(cor__SGZ_A1df))

hist(cor__SGZ_A1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A1 reactive Astrocyte gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_A1df <- cor__SGZ_A1df[order(cor__SGZ_A1df$correlation, decreasing = TRUE),]

fn_out_18 <- file.path(dir_outputs, "A1_reactive_Astrocyte_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_A1df,fn_out_18, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A2 astrocyte markers
#Calculat Spearman rank corr
cor_SGZ_A2 <- apply(logcounts(spe_SGZ_A2), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_A2 <- sapply(cor_SGZ_A2, function(x) x$estimate)
p_values_SGZ_A2 <- sapply(cor_SGZ_A2, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_A2df <- data.frame(gene_id = names(cor_SGZ_A2), correlation = cor_coeffs_SGZ_A2, adj.p.value = p_values_SGZ_A2)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_A2df) <- sub("\\.rho$", "", rownames(cor__SGZ_A2df))

hist(cor__SGZ_A2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A2 reactive Astrocyte gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_A2df <- cor__SGZ_A2df[order(cor__SGZ_A2df$correlation, decreasing = TRUE),]

fn_out_19 <- file.path(dir_outputs, "A2_reactive_Astrocyte_cor__SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_A2df,fn_out_19, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Ast1&6 astrocyte markers
#Calculat Spearman rank corr
cor_SGZ_Su <- apply(logcounts(spe_SGZ_Su), 1, function(row) cor.test(as.numeric(row), ages_SGZ, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_Su <- sapply(cor_SGZ_Su, function(x) x$estimate)
p_values_SGZ_Su <- sapply(cor_SGZ_Su, function(x) x$p.value)

# combine the results into a data frame
cor__SGZ_Sudf <- data.frame(gene_id = names(cor_SGZ_Su), correlation = cor_coeffs_SGZ_Su, adj.p.value = p_values_SGZ_Su)

# remove the ".rho" characters from the row names
rownames(cor__SGZ_Sudf) <- sub("\\.rho$", "", rownames(cor__SGZ_Sudf))

hist(cor__SGZ_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Ast1&6 reactive Astrocyte gene counts in SGZ")

# order the rows based on the correlation values
cor__SGZ_Sudf <- cor__SGZ_Sudf[order(cor__SGZ_Sudf$correlation, decreasing = TRUE),]

fn_out_20 <- file.path(dir_outputs, "Ast1&6_reactive_Astrocyte_SGZ_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__SGZ_Sudf,fn_out_20, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#PAN astrocyte markers
#Calculat Spearman rank corr
cor_GCL_PAN <- apply(logcounts(spe_GCL_PAN), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_PAN <- sapply(cor_GCL_PAN, function(x) x$estimate)
p_values_GCL_PAN <- sapply(cor_GCL_PAN, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_PANdf <- data.frame(gene_id = names(cor_GCL_PAN), correlation = cor_coeffs_GCL_PAN, adj.p.value = p_values_GCL_PAN)

# remove the ".rho" characters from the row names
rownames(cor__GCL_PANdf) <- sub("\\.rho$", "", rownames(cor__GCL_PANdf))

hist(cor__GCL_PANdf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "PAN reactive Astrocyte gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_PANdf <- cor__GCL_PANdf[order(cor__GCL_PANdf$correlation, decreasing = TRUE),]

fn_out_21 <- file.path(dir_outputs, "PAN_reactive_Astrocyte_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_PANdf,fn_out_21, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A1 astrocyte markers
#Calculat Spearman rank corr
cor_GCL_A1 <- apply(logcounts(spe_GCL_A1), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_A1 <- sapply(cor_GCL_A1, function(x) x$estimate)
p_values_GCL_A1 <- sapply(cor_GCL_A1, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_A1df <- data.frame(gene_id = names(cor_GCL_A1), correlation = cor_coeffs_GCL_A1, adj.p.value = p_values_GCL_A1)

# remove the ".rho" characters from the row names
rownames(cor__GCL_A1df) <- sub("\\.rho$", "", rownames(cor__GCL_A1df))

hist(cor__GCL_A1df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A1 reactive Astrocyte gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_A1df <- cor__GCL_A1df[order(cor__GCL_A1df$correlation, decreasing = TRUE),]

fn_out_22 <- file.path(dir_outputs, "A1_reactive_Astrocyte_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_A1df,fn_out_22, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#A2 astrocyte markers
#Calculat Spearman rank corr
cor_GCL_A2 <- apply(logcounts(spe_GCL_A2), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_GCL_A2 <- sapply(cor_GCL_A2, function(x) x$estimate)
p_values_GCL_A2 <- sapply(cor_GCL_A2, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_A2df <- data.frame(gene_id = names(cor_GCL_A2), correlation = cor_coeffs_GCL_A2, adj.p.value = p_values_GCL_A2)

# remove the ".rho" characters from the row names
rownames(cor__GCL_A2df) <- sub("\\.rho$", "", rownames(cor__GCL_A2df))

hist(cor__GCL_A2df$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "A2 reactive Astrocyte gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_A2df <- cor__GCL_A2df[order(cor__GCL_A2df$correlation, decreasing = TRUE),]

fn_out_23 <- file.path(dir_outputs, "A2_reactive_Astrocyte_cor__GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_A2df,fn_out_23, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
#Ast1&6 astrocyte markers
#Calculat Spearman rank corr
cor_GCL_Su <- apply(logcounts(spe_GCL_Su), 1, function(row) cor.test(as.numeric(row), ages_GCL, method = "spearman", adjust = "fdr"))

# extract the correlation coefficients and p-values from the results
cor_coeffs_SGZ_Su <- sapply(cor_GCL_Su, function(x) x$estimate)
p_values_SGZ_Su <- sapply(cor_GCL_Su, function(x) x$p.value)

# combine the results into a data frame
cor__GCL_Sudf <- data.frame(gene_id = names(cor_GCL_Su), correlation = cor_coeffs_SGZ_Su, adj.p.value = cor_coeffs_SGZ_Su)

# remove the ".rho" characters from the row names
rownames(cor__GCL_Sudf) <- sub("\\.rho$", "", rownames(cor__GCL_Sudf))

hist(cor__GCL_Sudf$correlation, xlab = "age rho correlation", ylab = "Gene count", main = "Ast1&6 reactive Astrocyte gene counts in GCL")

# order the rows based on the correlation values
cor__GCL_Sudf <- cor__GCL_Sudf[order(cor__GCL_Sudf$correlation, decreasing = TRUE),]

fn_out_24 <- file.path(dir_outputs, "Ast1&6_reactive_Astrocyte_GCL_rho_age_test_results")

# Export summary as .csv file
write.csv(cor__GCL_Sudf,fn_out_24, row.names = TRUE)
#---------------------------------------------------------------------------------------------------------------------------------------
