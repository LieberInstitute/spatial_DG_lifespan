#############################################
# spatial_DG_lifespan project
# Scatterpie plots with cell2location results
# Anthony Ramnauth, June 15 2023
#############################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(SpatialExperiment)
    library(dplyr)
    library(scatterpie)
    library(SPOTlight)
    library(spatialLIBD)
    library(viridis)
    library(ggcorrplot)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

#############################################################
# Scatterpie plots for deconvoluted cell types of Visium data
#############################################################

# Convert mean abundances to proportions

convert_to_proportions <- function(row) {
  proportions <- row / sum(row)
  return(proportions)
}

Br1412 <- spe[, spe$sample_id == "Br1412"]
Br1412y <- as.matrix(colData(Br1412)[, c(44:68)])
# Apply the function to each row of the data frame
Br1412y <- apply(Br1412y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br1412y <- as.data.frame(t(Br1412y))
for (col in 1:ncol(Br1412y)){
    colnames(Br1412y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br1412y)[col])
}

Br2706 <- spe[, spe$sample_id == "Br2706"]
Br2706y <- as.matrix(colData(Br2706)[, c(44:68)])
# Apply the function to each row of the data frame
Br2706y <- apply(Br2706y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br2706y <- as.data.frame(t(Br2706y))
for (col in 1:ncol(Br2706y)){
    colnames(Br2706y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br2706y)[col])
}

Br2720 <- spe[, spe$sample_id == "Br2720"]
Br2720y <- as.matrix(colData(Br2720)[, c(44:68)])
# Apply the function to each row of the data frame
Br2720y <- apply(Br2720y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br2720y <- as.data.frame(t(Br2720y))
for (col in 1:ncol(Br2720y)){
    colnames(Br2720y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br2720y)[col])
}

Br3942 <- spe[, spe$sample_id == "Br3942"]
Br3942y <- as.matrix(colData(Br3942)[, c(44:68)])
# Apply the function to each row of the data frame
Br3942y <- apply(Br3942y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br3942y <- as.data.frame(t(Br3942y))
for (col in 1:ncol(Br3942y)){
    colnames(Br3942y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br3942y)[col])
}

Br5242 <- spe[, spe$sample_id == "Br5242"]
Br5242y <- as.matrix(colData(Br5242)[, c(44:68)])
# Apply the function to each row of the data frame
Br5242y <- apply(Br5242y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br5242y <- as.data.frame(t(Br5242y))
for (col in 1:ncol(Br5242y)){
    colnames(Br5242y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br5242y)[col])
}

Br5699_new <- spe[, spe$sample_id == "Br5699_new"]
Br5699_newy <- as.matrix(colData(Br5699_new)[, c(44:68)])
# Apply the function to each row of the data frame
Br5699_newy <- apply(Br5699_newy, 1, convert_to_proportions)
# Convert the result back to a data frame
Br5699_newy <- as.data.frame(t(Br5699_newy))
for (col in 1:ncol(Br5699_newy)){
    colnames(Br5699_newy)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br5699_newy)[col])
}

Br6023 <- spe[, spe$sample_id == "Br6023"]
Br6023y <- as.matrix(colData(Br6023)[, c(44:68)])
# Apply the function to each row of the data frame
Br6023y <- apply(Br6023y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br6023y <- as.data.frame(t(Br6023y))
for (col in 1:ncol(Br6023y)){
    colnames(Br6023y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br6023y)[col])
}

Br6129_new <- spe[, spe$sample_id == "Br6129_new"]
Br6129_newy <- as.matrix(colData(Br6129_new)[, c(44:68)])
# Apply the function to each row of the data frame
Br6129_newy <- apply(Br6129_newy, 1, convert_to_proportions)
# Convert the result back to a data frame
Br6129_newy <- as.data.frame(t(Br6129_newy))
for (col in 1:ncol(Br6129_newy)){
    colnames(Br6129_newy)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br6129_newy)[col])
}

Br6299_new <- spe[, spe$sample_id == "Br6299_new"]
Br6299_newy <- as.matrix(colData(Br6299_new)[, c(44:68)])
# Apply the function to each row of the data frame
Br6299_newy <- apply(Br6299_newy, 1, convert_to_proportions)
# Convert the result back to a data frame
Br6299_newy <- as.data.frame(t(Br6299_newy))
for (col in 1:ncol(Br6299_newy)){
    colnames(Br6299_newy)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br6299_newy)[col])
}

Br6522 <- spe[, spe$sample_id == "Br6522"]
Br6522y <- as.matrix(colData(Br6522)[, c(44:68)])
# Apply the function to each row of the data frame
Br6522y <- apply(Br6522y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br6522y <- as.data.frame(t(Br6522y))
for (col in 1:ncol(Br6522y)){
    colnames(Br6522y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br6522y)[col])
}

Br8181 <- spe[, spe$sample_id == "Br8181"]
Br8181y <- as.matrix(colData(Br8181)[, c(44:68)])
# Apply the function to each row of the data frame
Br8181y <- apply(Br8181y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br8181y <- as.data.frame(t(Br8181y))
for (col in 1:ncol(Br8181y)){
    colnames(Br8181y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br8181y)[col])
}

Br8195 <- spe[, spe$sample_id == "Br8195"]
Br8195y <- as.matrix(colData(Br8195)[, c(44:68)])
# Apply the function to each row of the data frame
Br8195y <- apply(Br8195y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br8195y <- as.data.frame(t(Br8195y))
for (col in 1:ncol(Br8195y)){
    colnames(Br8195y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br8195y)[col])
}

Br8533 <- spe[, spe$sample_id == "Br8533"]
Br8533y <- as.matrix(colData(Br8533)[, c(44:68)])
# Apply the function to each row of the data frame
Br8533y <- apply(Br8533y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br8533y <- as.data.frame(t(Br8533y))
for (col in 1:ncol(Br8533y)){
    colnames(Br8533y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br8533y)[col])
}

Br8667 <- spe[, spe$sample_id == "Br8667"]
Br8667y <- as.matrix(colData(Br8667)[, c(44:68)])
# Apply the function to each row of the data frame
Br8667y <- apply(Br8667y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br8667y <- as.data.frame(t(Br8667y))
for (col in 1:ncol(Br8667y)){
    colnames(Br8667y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br8667y)[col])
}

Br8686 <- spe[, spe$sample_id == "Br8686"]
Br8686y <- as.matrix(colData(Br8686)[, c(44:68)])
# Apply the function to each row of the data frame
Br8686y <- apply(Br8686y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br8686y <- as.data.frame(t(Br8686y))
for (col in 1:ncol(Br8686y)){
    colnames(Br8686y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br8686y)[col])
}

Br8700 <- spe[, spe$sample_id == "Br8700"]
Br8700y <- as.matrix(colData(Br8700)[, c(44:68)])
# Apply the function to each row of the data frame
Br8700y <- apply(Br8700y, 1, convert_to_proportions)
# Convert the result back to a data frame
Br8700y <- as.data.frame(t(Br8700y))
for (col in 1:ncol(Br8700y)){
    colnames(Br8700y)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(Br8700y)[col])
}

cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3",
    "VLMC" = "firebrick1", "Pericyte" = "red2", "Endoth" = "red", "SMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "yellow", "Astro_2" = "yellow3")


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell2loc_scatterpie_plots.pdf"), width = 12, height = 8)

plotSpatialScatterpie(
    x = Br1412,
    y = Br1412y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br1412 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br2706,
    y = Br2706y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br2706 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br2720,
    y = Br2720y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br2720 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br3942,
    y = Br3942y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br3942 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br5242,
    y = Br5242y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br5242 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br5699_new,
    y = Br5699_newy,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br5699_new cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br6023,
    y = Br6023y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br6023 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br6129_new,
    y = Br6129_newy,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br6129_new cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br6299_new,
    y = Br6299_newy,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br6299_new cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br6522,
    y = Br6522y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br6522 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8181,
    y = Br8181y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8181 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8195,
    y = Br8195y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8195 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8533,
    y = Br8533y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8533 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8667,
    y = Br8667y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8667 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8686,
    y = Br8686y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8686 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8700,
    y = Br8700y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8700 cell-type deconvolution") +
    scale_y_reverse()

dev.off()

