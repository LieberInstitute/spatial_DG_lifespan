###########################################
# spatial_DG_lifespan project
# MDS plot of Cell-Type spatial similarity
# Anthony Ramnauth, Oct 20 2022
###########################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(SpatialExperiment)
    library(dplyr)
    library(CARD)
    library(viridis)
    library(spatialLIBD)
    library(smacof)
    library(ggplot2)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

############################################################################
# Make correlation matrix for cell-type proportions across spatial locations
############################################################################

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

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

# All Visium spots

spe_cells <- as.data.frame(colData(spe)[, c(68:97)])

#reorder correlation by column index

####################################################
# Run Multi-dimensional scaling for all Visium spots
####################################################

# All
Rmat <- cor(spe_cells, use = "pairwise.complete.obs")

Delta <- sim2diss(Rmat, method = "corr", to.dist = TRUE)

round(Delta, 2)

mds_1 <- mds(Delta)

mds_2 <- mds(Delta, type = "ordinal")
plot(mds_1, plot.type = "Shepard", main = "Shepard Diagram (Ratio Transformation)")

plot(mds_2, main = "Spatial Similarity of Cell-types")

cell_sim <- as.data.frame(mds_2$conf)
cell_sim$type <- factor(rownames(cell_sim), levels = unique(rownames(cell_sim)))

# Age groups for all Visium spots
#infant
infant_spe <- spe[, spe$age_bin == "Infant"]
infant_cells <- as.data.frame(colData(infant_spe)[, c(44:66)])
infant_cells <- infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_infant_cells = cor(infant_cells, use = "pairwise.complete.obs")
infant_Delta <- sim2diss(cor_infant_cells, method = "corr", to.dist = TRUE)
infant_mds_1 <- mds(infant_Delta)
infant_mds_2 <- mds(infant_Delta, type = "ordinal")
infant_cell_sim <- as.data.frame(infant_mds_2$conf)
infant_cell_sim$type <- factor(rownames(infant_cell_sim), levels = unique(rownames(infant_cell_sim)))

#teen
teen_spe <- spe[, spe$age_bin == "Teen"]
teen_cells <- as.data.frame(colData(teen_spe)[, c(44:66)])
teen_cells <- teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_teen_cells = cor(teen_cells, use = "pairwise.complete.obs")
teen_Delta <- sim2diss(cor_teen_cells, method = "corr", to.dist = TRUE)
teen_mds_1 <- mds(teen_Delta)
teen_mds_2 <- mds(teen_Delta, type = "ordinal")
teen_cell_sim <- as.data.frame(teen_mds_2$conf)
teen_cell_sim$type <- factor(rownames(teen_cell_sim), levels = unique(rownames(teen_cell_sim)))

#adult
adult_spe <- spe[, spe$age_bin == "Adult"]
adult_cells <- as.data.frame(colData(adult_spe)[, c(44:66)])
adult_cells <- adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_adult_cells = cor(adult_cells, use = "pairwise.complete.obs")
adult_Delta <- sim2diss(cor_adult_cells, method = "corr", to.dist = TRUE)
adult_mds_1 <- mds(adult_Delta)
adult_mds_2 <- mds(adult_Delta, type = "ordinal")
adult_cell_sim <- as.data.frame(adult_mds_2$conf)
adult_cell_sim$type <- factor(rownames(adult_cell_sim), levels = unique(rownames(adult_cell_sim)))

#elderly
elderly_spe <- spe[, spe$age_bin == "Elderly"]
elderly_cells <- as.data.frame(colData(elderly_spe)[, c(44:66)])
elderly_cells <- elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_elderly_cells = cor(elderly_cells, use = "pairwise.complete.obs")
elderly_Delta <- sim2diss(cor_elderly_cells, method = "corr", to.dist = TRUE)
elderly_mds_1 <- mds(elderly_Delta)
elderly_mds_2 <- mds(elderly_Delta, type = "ordinal")
elderly_cell_sim <- as.data.frame(elderly_mds_2$conf)
elderly_cell_sim$type <- factor(rownames(elderly_cell_sim), levels = unique(rownames(elderly_cell_sim)))


###########################
# Plot for All Visium Spots
###########################

cell_colors <- c("Oligo_1" = "plum3", "Oligo_2" = "plum4", "Microglia" = "tan2", "Macrocyte" = "tan1",
    "OPC_1" = "goldenrod", "OPC_2" = "goldenrod3", "InN_LAMP5" = "green", "InN_VIP" = "green1",
    "InN_SST" = "green2", "InN_PV" = "green3", "InN_NR2F2" = "green4", "InN_LHX6" = "lawngreen",
    "InN_MEIS2" = "mediumseagreen", "Cajal_Ret" = "black", "Vasc_LM" = "red", "Artl_S_Muscle" = "red1",
    "Pericyte" = "red2", "Endoth" = "red3", "Vasc_S_Muscle" = "red4", "T_cell" = "tan3",
    "Myeloid" = "tan4", "COP" = "goldenrod4", "GC" = "blue", "CA3_N" = "dodgerblue", "EC_N" = "blue1",
    "Mossy" = "blue4", "CA1_N" = "blue2", "SUB_N" = "blue3", "Astro_1" = "yellow2", "Astro_2" = "yellow")

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_MDS_spots.pdf"), width = 12,
    height = 8)

ggplot(cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(infant_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(infant_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Infant Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(teen_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(teen_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Teen Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(adult_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(adult_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Adult Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(elderly_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(elderly_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Elderly Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

dev.off()

#######################################################
# Run Multi-dimensional scaling for BayesSpace clusters
#######################################################

#Cluster 1
spe_cluster1 <- spe[, spe$bayesSpace_harmony_8 == 1]
cluster1_cells <- as.data.frame(colData(spe_cluster1)[, c(44:66)])
cluster1_cells <- cluster1_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_cells = cor(cluster1_cells, use = "pairwise.complete.obs")
cluster1_Delta <- sim2diss(cor_cluster1_cells, method = "corr", to.dist = TRUE)
cluster1_mds_1 <- mds(cluster1_Delta)
cluster1_mds_2 <- mds(cluster1_Delta, type = "ordinal")
cluster1_cell_sim <- as.data.frame(cluster1_mds_2$conf)
cluster1_cell_sim$type <- factor(rownames(cluster1_cell_sim), levels = unique(rownames(cluster1_cell_sim)))

#infant
infant_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Infant"]
cluster1_infant_cells <- as.data.frame(colData(infant_cluster1)[, c(44:66)])
cluster1_infant_cells <- cluster1_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_infant_cells = cor(cluster1_infant_cells, use = "pairwise.complete.obs")
cluster1_infant_Delta <- sim2diss(cor_cluster1_infant_cells, method = "corr", to.dist = TRUE)
cluster1_infant_mds_1 <- mds(cluster1_infant_Delta)
cluster1_infant_mds_2 <- mds(cluster1_infant_Delta, type = "ordinal")
cluster1_infant_cell_sim <- as.data.frame(cluster1_infant_mds_2$conf)
cluster1_infant_cell_sim$type <- factor(rownames(cluster1_infant_cell_sim),
    levels = unique(rownames(cluster1_infant_cell_sim)))

#teen
teen_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Teen"]
cluster1_teen_cells <- as.data.frame(colData(teen_cluster1)[, c(44:66)])
cluster1_teen_cells <- cluster1_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_teen_cells = cor(cluster1_teen_cells, use = "pairwise.complete.obs")
cluster1_teen_Delta <- sim2diss(cor_cluster1_teen_cells, method = "corr", to.dist = TRUE)
cluster1_teen_mds_1 <- mds(cluster1_teen_Delta)
cluster1_teen_mds_2 <- mds(cluster1_teen_Delta, type = "ordinal")
cluster1_teen_cell_sim <- as.data.frame(cluster1_teen_mds_2$conf)
cluster1_teen_cell_sim$type <- factor(rownames(cluster1_teen_cell_sim),
    levels = unique(rownames(cluster1_teen_cell_sim)))

#adult
adult_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Adult"]
cluster1_adult_cells <- as.data.frame(colData(adult_cluster1)[, c(44:66)])
cluster1_adult_cells <- cluster1_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_adult_cells = cor(cluster1_adult_cells, use = "pairwise.complete.obs")
cluster1_adult_Delta <- sim2diss(cor_cluster1_adult_cells, method = "corr", to.dist = TRUE)
cluster1_adult_mds_1 <- mds(cluster1_adult_Delta)
cluster1_adult_mds_2 <- mds(cluster1_adult_Delta, type = "ordinal")
cluster1_adult_cell_sim <- as.data.frame(cluster1_adult_mds_2$conf)
cluster1_adult_cell_sim$type <- factor(rownames(cluster1_adult_cell_sim),
    levels = unique(rownames(cluster1_adult_cell_sim)))

#elderly
elderly_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Elderly"]
cluster1_elderly_cells <- as.data.frame(colData(elderly_cluster1)[, c(44:66)])
cluster1_elderly_cells <- cluster1_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_elderly_cells = cor(cluster1_elderly_cells, use = "pairwise.complete.obs")
cluster1_elderly_Delta <- sim2diss(cor_cluster1_elderly_cells, method = "corr", to.dist = TRUE)
cluster1_elderly_mds_1 <- mds(cluster1_elderly_Delta)
cluster1_elderly_mds_2 <- mds(cluster1_elderly_Delta, type = "ordinal")
cluster1_elderly_cell_sim <- as.data.frame(cluster1_elderly_mds_2$conf)
cluster1_elderly_cell_sim$type <- factor(rownames(cluster1_elderly_cell_sim),
    levels = unique(rownames(cluster1_elderly_cell_sim)))

#Cluster 2
spe_cluster2 <- spe[, spe$bayesSpace_harmony_8 == 2]
cluster2_cells <- as.data.frame(colData(spe_cluster2)[, c(44:66)])
cluster2_cells <- cluster2_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_cells = cor(cluster2_cells, use = "pairwise.complete.obs")
cluster2_Delta <- sim2diss(cor_cluster2_cells, method = "corr", to.dist = TRUE)
cluster2_mds_1 <- mds(cluster2_Delta)
cluster2_mds_2 <- mds(cluster2_Delta, type = "ordinal")
cluster2_cell_sim <- as.data.frame(cluster2_mds_2$conf)
cluster2_cell_sim$type <- factor(rownames(cluster2_cell_sim), levels = unique(rownames(cluster2_cell_sim)))

#infant
infant_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Infant"]
cluster2_infant_cells <- as.data.frame(colData(infant_cluster2)[, c(44:66)])
cluster2_infant_cells <- cluster2_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_infant_cells = cor(cluster2_infant_cells, use = "pairwise.complete.obs")
cluster2_infant_Delta <- sim2diss(cor_cluster2_infant_cells, method = "corr", to.dist = TRUE)
cluster2_infant_mds_1 <- mds(cluster2_infant_Delta)
cluster2_infant_mds_2 <- mds(cluster2_infant_Delta, type = "ordinal")
cluster2_infant_cell_sim <- as.data.frame(cluster2_infant_mds_2$conf)
cluster2_infant_cell_sim$type <- factor(rownames(cluster2_infant_cell_sim),
    levels = unique(rownames(cluster2_infant_cell_sim)))

#teen
teen_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Teen"]
cluster2_teen_cells <- as.data.frame(colData(teen_cluster2)[, c(44:66)])
cluster2_teen_cells <- cluster2_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_teen_cells = cor(cluster2_teen_cells, use = "pairwise.complete.obs")
cluster2_teen_Delta <- sim2diss(cor_cluster2_teen_cells, method = "corr", to.dist = TRUE)
cluster2_teen_mds_1 <- mds(cluster2_teen_Delta)
cluster2_teen_mds_2 <- mds(cluster2_teen_Delta, type = "ordinal")
cluster2_teen_cell_sim <- as.data.frame(cluster2_teen_mds_2$conf)
cluster2_teen_cell_sim$type <- factor(rownames(cluster2_teen_cell_sim),
    levels = unique(rownames(cluster2_teen_cell_sim)))

#adult
adult_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Adult"]
cluster2_adult_cells <- as.data.frame(colData(adult_cluster2)[, c(44:66)])
cluster2_adult_cells <- cluster2_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_adult_cells = cor(cluster2_adult_cells, use = "pairwise.complete.obs")
cluster2_adult_Delta <- sim2diss(cor_cluster2_adult_cells, method = "corr", to.dist = TRUE)
cluster2_adult_mds_1 <- mds(cluster2_adult_Delta)
cluster2_adult_mds_2 <- mds(cluster2_adult_Delta, type = "ordinal")
cluster2_adult_cell_sim <- as.data.frame(cluster2_adult_mds_2$conf)
cluster2_adult_cell_sim$type <- factor(rownames(cluster2_adult_cell_sim),
    levels = unique(rownames(cluster2_adult_cell_sim)))

#elderly
elderly_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Elderly"]
cluster2_elderly_cells <- as.data.frame(colData(elderly_cluster2)[, c(44:66)])
cluster2_elderly_cells <- cluster2_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_elderly_cells = cor(cluster2_elderly_cells, use = "pairwise.complete.obs")
cluster2_elderly_Delta <- sim2diss(cor_cluster2_elderly_cells, method = "corr", to.dist = TRUE)
cluster2_elderly_mds_1 <- mds(cluster2_elderly_Delta)
cluster2_elderly_mds_2 <- mds(cluster2_elderly_Delta, type = "ordinal")
cluster2_elderly_cell_sim <- as.data.frame(cluster2_elderly_mds_2$conf)
cluster2_elderly_cell_sim$type <- factor(rownames(cluster2_elderly_cell_sim),
    levels = unique(rownames(cluster2_elderly_cell_sim)))

#Cluster 4
spe_cluster4 <- spe[, spe$bayesSpace_harmony_8 == 4]
cluster4_cells <- as.data.frame(colData(spe_cluster4)[, c(44:66)])
cluster4_cells <- cluster4_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_cells = cor(cluster4_cells, use = "pairwise.complete.obs")
cluster4_Delta <- sim2diss(cor_cluster4_cells, method = "corr", to.dist = TRUE)
cluster4_mds_1 <- mds(cluster4_Delta)
cluster4_mds_2 <- mds(cluster4_Delta, type = "ordinal")
cluster4_cell_sim <- as.data.frame(cluster4_mds_2$conf)
cluster4_cell_sim$type <- factor(rownames(cluster4_cell_sim), levels = unique(rownames(cluster4_cell_sim)))

#infant
infant_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Infant"]
cluster4_infant_cells <- as.data.frame(colData(infant_cluster4)[, c(44:66)])
cluster4_infant_cells <- cluster4_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_infant_cells = cor(cluster4_infant_cells, use = "pairwise.complete.obs")
cluster4_infant_Delta <- sim2diss(cor_cluster4_infant_cells, method = "corr", to.dist = TRUE)
cluster4_infant_mds_1 <- mds(cluster4_infant_Delta)
cluster4_infant_mds_2 <- mds(cluster4_infant_Delta, type = "ordinal")
cluster4_infant_cell_sim <- as.data.frame(cluster4_infant_mds_2$conf)
cluster4_infant_cell_sim$type <- factor(rownames(cluster4_infant_cell_sim),
    levels = unique(rownames(cluster4_infant_cell_sim)))

#teen
teen_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Teen"]
cluster4_teen_cells <- as.data.frame(colData(teen_cluster4)[, c(44:66)])
cluster4_teen_cells <- cluster4_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_teen_cells = cor(cluster4_teen_cells, use = "pairwise.complete.obs")
cluster4_teen_Delta <- sim2diss(cor_cluster4_teen_cells, method = "corr", to.dist = TRUE)
cluster4_teen_mds_1 <- mds(cluster4_teen_Delta)
cluster4_teen_mds_2 <- mds(cluster4_teen_Delta, type = "ordinal")
cluster4_teen_cell_sim <- as.data.frame(cluster4_teen_mds_2$conf)
cluster4_teen_cell_sim$type <- factor(rownames(cluster4_teen_cell_sim),
    levels = unique(rownames(cluster4_teen_cell_sim)))

#adult
adult_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Adult"]
cluster4_adult_cells <- as.data.frame(colData(adult_cluster4)[, c(44:66)])
cluster4_adult_cells <- cluster4_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_adult_cells = cor(cluster4_adult_cells, use = "pairwise.complete.obs")
cluster4_adult_Delta <- sim2diss(cor_cluster4_adult_cells, method = "corr", to.dist = TRUE)
cluster4_adult_mds_1 <- mds(cluster4_adult_Delta)
cluster4_adult_mds_2 <- mds(cluster4_adult_Delta, type = "ordinal")
cluster4_adult_cell_sim <- as.data.frame(cluster4_adult_mds_2$conf)
cluster4_adult_cell_sim$type <- factor(rownames(cluster4_adult_cell_sim),
    levels = unique(rownames(cluster4_adult_cell_sim)))

#elderly
elderly_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Elderly"]
cluster4_elderly_cells <- as.data.frame(colData(elderly_cluster4)[, c(44:66)])
cluster4_elderly_cells <- cluster4_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_elderly_cells = cor(cluster4_elderly_cells, use = "pairwise.complete.obs")
cluster4_elderly_Delta <- sim2diss(cor_cluster4_elderly_cells, method = "corr", to.dist = TRUE)
cluster4_elderly_mds_1 <- mds(cluster4_elderly_Delta)
cluster4_elderly_mds_2 <- mds(cluster4_elderly_Delta, type = "ordinal")
cluster4_elderly_cell_sim <- as.data.frame(cluster4_elderly_mds_2$conf)
cluster4_elderly_cell_sim$type <- factor(rownames(cluster4_elderly_cell_sim),
    levels = unique(rownames(cluster4_elderly_cell_sim)))

#Cluster 8
spe_cluster8 <- spe[, spe$bayesSpace_harmony_8 == 8]
cluster8_cells <- as.data.frame(colData(spe_cluster8)[, c(44:66)])
cluster8_cells <- cluster8_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_cells = cor(cluster8_cells, use = "pairwise.complete.obs")
cluster8_Delta <- sim2diss(cor_cluster8_cells, method = "corr", to.dist = TRUE)
cluster8_mds_1 <- mds(cluster8_Delta)
cluster8_mds_2 <- mds(cluster8_Delta, type = "ordinal")
cluster8_cell_sim <- as.data.frame(cluster8_mds_2$conf)
cluster8_cell_sim$type <- factor(rownames(cluster8_cell_sim), levels = unique(rownames(cluster8_cell_sim)))

#infant
infant_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Infant"]
cluster8_infant_cells <- as.data.frame(colData(infant_cluster8)[, c(44:66)])
cluster8_infant_cells <- cluster8_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_infant_cells = cor(cluster8_infant_cells, use = "pairwise.complete.obs")
cluster8_infant_Delta <- sim2diss(cor_cluster8_infant_cells, method = "corr", to.dist = TRUE)
cluster8_infant_mds_1 <- mds(cluster8_infant_Delta)
cluster8_infant_mds_2 <- mds(cluster8_infant_Delta, type = "ordinal")
cluster8_infant_cell_sim <- as.data.frame(cluster8_infant_mds_2$conf)
cluster8_infant_cell_sim$type <- factor(rownames(cluster8_infant_cell_sim),
    levels = unique(rownames(cluster8_infant_cell_sim)))

#teen
teen_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Teen"]
cluster8_teen_cells <- as.data.frame(colData(teen_cluster8)[, c(44:66)])
cluster8_teen_cells <- cluster8_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_teen_cells = cor(cluster8_teen_cells, use = "pairwise.complete.obs")
cluster8_teen_Delta <- sim2diss(cor_cluster8_teen_cells, method = "corr", to.dist = TRUE)
cluster8_teen_mds_1 <- mds(cluster8_teen_Delta)
cluster8_teen_mds_2 <- mds(cluster8_teen_Delta, type = "ordinal")
cluster8_teen_cell_sim <- as.data.frame(cluster8_teen_mds_2$conf)
cluster8_teen_cell_sim$type <- factor(rownames(cluster8_teen_cell_sim),
    levels = unique(rownames(cluster8_teen_cell_sim)))

#adult
adult_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Adult"]
cluster8_adult_cells <- as.data.frame(colData(adult_cluster8)[, c(44:66)])
cluster8_adult_cells <- cluster8_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_adult_cells = cor(cluster8_adult_cells, use = "pairwise.complete.obs")
cluster8_adult_Delta <- sim2diss(cor_cluster8_adult_cells, method = "corr", to.dist = TRUE)
cluster8_adult_mds_1 <- mds(cluster8_adult_Delta)
cluster8_adult_mds_2 <- mds(cluster8_adult_Delta, type = "ordinal")
cluster8_adult_cell_sim <- as.data.frame(cluster8_adult_mds_2$conf)
cluster8_adult_cell_sim$type <- factor(rownames(cluster8_adult_cell_sim),
    levels = unique(rownames(cluster8_adult_cell_sim)))

#elderly
elderly_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Elderly"]
cluster8_elderly_cells <- as.data.frame(colData(elderly_cluster8)[, c(44:66)])
cluster8_elderly_cells <- cluster8_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_elderly_cells = cor(cluster8_elderly_cells, use = "pairwise.complete.obs")
cluster8_elderly_Delta <- sim2diss(cor_cluster8_elderly_cells, method = "corr", to.dist = TRUE)
cluster8_elderly_mds_1 <- mds(cluster8_elderly_Delta)
cluster8_elderly_mds_2 <- mds(cluster8_elderly_Delta, type = "ordinal")
cluster8_elderly_cell_sim <- as.data.frame(cluster8_elderly_mds_2$conf)
cluster8_elderly_cell_sim$type <- factor(rownames(cluster8_elderly_cell_sim),
    levels = unique(rownames(cluster8_elderly_cell_sim)))

###################################
# MDS Plots for BayesSpace clusters
###################################

# Plot for cluster 1

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_MDS_cluster1.pdf"), width = 12,
    height = 8)

ggplot(cluster1_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster1_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_1 spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster1_infant_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster1_infant_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_1 Infant Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster1_teen_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster1_teen_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_1 Teen Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster1_adult_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster1_adult_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_1 Adult Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster1_elderly_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster1_elderly_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types all BayesSpace cluster_1 across Elderly Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

dev.off()

# Plot for cluster 2

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_MDS_cluster2.pdf"), width = 12,
    height = 8)

ggplot(cluster2_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster2_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_2 spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster2_infant_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster2_infant_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_2 Infant Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster2_teen_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster2_teen_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_2 Teen Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster2_adult_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster2_adult_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_2 Adult Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster2_elderly_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster2_elderly_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types all BayesSpace cluster_2 across Elderly Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

dev.off()

# Plot for cluster 4

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_MDS_cluster4.pdf"), width = 12,
    height = 8)

ggplot(cluster4_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster4_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_4 spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster4_infant_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster4_infant_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_4 Infant Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster4_teen_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster4_teen_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_4 Teen Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster4_adult_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster4_adult_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_4 Adult Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster4_elderly_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster4_elderly_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types all BayesSpace cluster_4 across Elderly Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

dev.off()

# Plot for cluster 8

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_MDS_cluster8.pdf"), width = 12,
    height = 8)

ggplot(cluster8_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster8_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_8 spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster8_infant_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster8_infant_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_8 Infant Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster8_teen_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster8_teen_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_8 Teen Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster8_adult_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster8_adult_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all BayesSpace cluster_8 Adult Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(cluster8_elderly_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cluster8_elderly_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    coord_fixed() +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types all BayesSpace cluster_8 across Elderly Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
