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
    library(viridis)
    library(spatialLIBD)
    library(smacof)
    library(ggplot2)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

############################################################################
# Make correlation matrix for cell-type proportions across spatial locations
############################################################################

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        between(spe.age, 0, 3) ~ "Infant",
        between(spe.age, 13, 19) ~ "Teen",
        between(spe.age, 20, 50) ~ "Adult",
        between(spe.age, 60, 100) ~ "Elderly"
    ))

stopifnot(age_df$spe.key == spe$key)

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

# All Visium spots

spe_cells <- as.data.frame(colData(spe)[, c(44:68)])

for ( col in 1:ncol(spe_cells)){
    colnames(spe_cells)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(spe_cells)[col])
}

# reorder correlation by column index

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
plot(mds_2, plot.type = "Shepard", main = "Shepard Diagram (Ordinal Transformation)")
plot(mds_2, main = "Spatial Similarity of Cell-types")

cell_sim <- as.data.frame(mds_2$conf)
cell_sim$type <- factor(rownames(cell_sim), levels = unique(rownames(cell_sim)))

# Age groups for all Visium spots
#infant
infant_spe <- spe[, spe$age_bin == "Infant"]
infant_cells <- as.data.frame(colData(infant_spe)[, c(44:68)])
for ( col in 1:ncol(infant_cells)){
    colnames(infant_cells)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(infant_cells)[col])
}
cor_infant_cells = cor(infant_cells, use = "pairwise.complete.obs")
infant_Delta <- sim2diss(cor_infant_cells, method = "corr", to.dist = TRUE)
infant_mds_1 <- mds(infant_Delta)
infant_mds_2 <- mds(infant_Delta, type = "ordinal")
infant_cell_sim <- as.data.frame(infant_mds_2$conf)
infant_cell_sim$type <- factor(rownames(infant_cell_sim), levels = unique(rownames(infant_cell_sim)))

#teen
teen_spe <- spe[, spe$age_bin == "Teen"]
teen_cells <- as.data.frame(colData(teen_spe)[, c(44:68)])
for ( col in 1:ncol(teen_cells)){
    colnames(teen_cells)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(teen_cells)[col])
}
cor_teen_cells = cor(teen_cells, use = "pairwise.complete.obs")
teen_Delta <- sim2diss(cor_teen_cells, method = "corr", to.dist = TRUE)
teen_mds_1 <- mds(teen_Delta)
teen_mds_2 <- mds(teen_Delta, type = "ordinal")
teen_cell_sim <- as.data.frame(teen_mds_2$conf)
teen_cell_sim$type <- factor(rownames(teen_cell_sim), levels = unique(rownames(teen_cell_sim)))

#adult
adult_spe <- spe[, spe$age_bin == "Adult"]
adult_cells <- as.data.frame(colData(adult_spe)[, c(44:68)])
for ( col in 1:ncol(adult_cells)){
    colnames(adult_cells)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(adult_cells)[col])
}
cor_adult_cells = cor(adult_cells, use = "pairwise.complete.obs")
adult_Delta <- sim2diss(cor_adult_cells, method = "corr", to.dist = TRUE)
adult_mds_1 <- mds(adult_Delta)
adult_mds_2 <- mds(adult_Delta, type = "ordinal")
adult_cell_sim <- as.data.frame(adult_mds_2$conf)
adult_cell_sim$type <- factor(rownames(adult_cell_sim), levels = unique(rownames(adult_cell_sim)))

#elderly
elderly_spe <- spe[, spe$age_bin == "Elderly"]
elderly_cells <- as.data.frame(colData(elderly_spe)[, c(44:68)])
for ( col in 1:ncol(elderly_cells)){
    colnames(elderly_cells)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(elderly_cells)[col])
}
cor_elderly_cells = cor(elderly_cells, use = "pairwise.complete.obs")
elderly_Delta <- sim2diss(cor_elderly_cells, method = "corr", to.dist = TRUE)
elderly_mds_1 <- mds(elderly_Delta)
elderly_mds_2 <- mds(elderly_Delta, type = "ordinal")
elderly_cell_sim <- as.data.frame(elderly_mds_2$conf)
elderly_cell_sim$type <- factor(rownames(elderly_cell_sim), levels = unique(rownames(elderly_cell_sim)))


###########################
# Plot for All Visium Spots
###########################

cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3",
    "VLMC" = "firebrick1", "Pericyte" = "red2", "Endoth" = "red", "SMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "yellow", "Astro_2" = "yellow3")

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_MDS_spots.pdf"), width = 12,
    height = 8)

ggplot(cell_sim, aes(x = D1, y = D2, color = type, label = rownames(cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across all Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(infant_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(infant_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Infant Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(teen_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(teen_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Teen Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(adult_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(adult_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Adult Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

ggplot(elderly_cell_sim, aes(x = D1, y = D2, color = type, label = rownames(elderly_cell_sim))) +
    geom_point(size = 5) +
    geom_text(size = 3.5, vjust = -0.8) +
    scale_colour_manual(values = cell_colors) +
    ggtitle("MDS Plot of Spatial Similarity for Cell-types across Elderly Visium spots") +
    labs(x = "MDS1", y = "MDS2") +
    theme_bw() +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"))

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
