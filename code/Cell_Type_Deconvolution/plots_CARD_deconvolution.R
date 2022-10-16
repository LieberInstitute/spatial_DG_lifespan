############################################
# spatial_DG_lifespan project
# Plotting Cell-Type Deconvolution with CARD
# Anthony Ramnauth, Oct 15 2022
############################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(SpatialExperiment)
    library(dplyr)
    library(scatterpie)
    library(SPOTlight)
    library(spatialLIBD)
    library(CARD)
    library(viridis)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

##################################################
# Plots for deconvoluted cell types of Visium data
##################################################

Br1412 <- spe[, spe$sample_id == "Br1412"]
Br1412y <- as.matrix(colData(Br1412)[, c(44:66)])

Br2706 <- spe[, spe$sample_id == "Br2706"]
Br2706y <- as.matrix(colData(Br2706)[, c(44:66)])

Br3942 <- spe[, spe$sample_id == "Br3942"]
Br3942y <- as.matrix(colData(Br3942)[, c(44:66)])

Br5242 <- spe[, spe$sample_id == "Br5242"]
Br5242y <- as.matrix(colData(Br5242)[, c(44:66)])

Br6023 <- spe[, spe$sample_id == "Br6023"]
Br6023y <- as.matrix(colData(Br6023)[, c(44:66)])

Br8195 <- spe[, spe$sample_id == "Br8195"]
Br8195y <- as.matrix(colData(Br8195)[, c(44:66)])

Br8667 <- spe[, spe$sample_id == "Br8667"]
Br8667y <- as.matrix(colData(Br8667)[, c(44:66)])

Br8686 <- spe[, spe$sample_id == "Br8686"]
Br8686y <- as.matrix(colData(Br8686)[, c(44:66)])

cell_colors <- c("Oligo" = "plum3", "Micro" = "tan2", "OPC" = "goldenrod", "Inhib_A" = "green",
    "Astro_A" = "yellow2", "Excit_B" = "dodgerblue", "Astro_B" = "yellow", "Excit_C" = "dodgerblue1",
    "Tcell"= "darkred", "drop.doublet" = "black", "Inhib_B" = "green2", "Excit_H" = "dodgerblue2",
    "drop.lowNTx_B" = "white", "drop.lowNTx_A" = "ghostwhite", "Excit_G" = "dodgerblue3",
    "Mural" = "brown", "Excit_F" = "dodgerblue4", "OPC_COP" = "goldenrod4", "Inhib_C" = "green3",
    "Excit_A" = "blue2", "Inhib_D" = "green4", "Excit_E" = "blue4",
    "Excit_D" = "midnightblue")

pdf(file = here::here("Github", "spatial_DG_lifespan", "plots", "Cell_Type_Deconvolution",
    "Scatterpie_plots.pdf"))
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
    x = Br6023,
    y = Br6023y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br6023 cell-type deconvolution") +
    scale_y_reverse()

plotSpatialScatterpie(
    x = Br8195,
    y = Br8195y,
) +
    scale_fill_manual(values = cell_colors) +
    labs(title = "Br8195 cell-type deconvolution") +
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

dev.off()

#############################################################
# Plot cell-type proportions separately for each capture area
#############################################################

spe_cells <- as.data.frame(cbind(colData(spe)$sample_id, colData(spe)[, c(44:66)], spatialCoords(spe)))
spe_cells <- rename(spe_cells, "sample_id" = "colData.spe..sample_id")

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Cell_type_full_scale_proportions.pdf"))

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Oligo)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Oligo_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Micro)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Micro_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = OPC)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "OPC_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_A)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Inhib_A_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Astro_A)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Astro_A_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_B)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_B_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Astro_B)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Astro_B_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_C)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_C_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Tcell)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Tcell_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = drop.doublet)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "drop.doublet__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_B)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Inhib_B_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_H)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_H_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = drop.lowNTx_B)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "drop.lowNTx_B__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = drop.lowNTx_A)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "drop.lowNTx_A__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_G)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_G_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Mural)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Mural_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_F)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_F_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = OPC_COP)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "OPC_COP_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_C)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Inhib_C_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_A)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_A_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_D)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Inhib_D_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_E)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_E_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_D)) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Excit_D_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

dev.off()

###################################################################################################
# Replotting with relative proportion scales for better highlighting of relative spatial expression
###################################################################################################

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Cell_type_relative_scale_proportions.pdf"))

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Oligo)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Oligo_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Micro)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Micro_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = OPC)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "OPC_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_A)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Inhib_A_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Astro_A)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Astro_A_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_B)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_B_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Astro_B)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Astro_B_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_C)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_C_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Tcell)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Tcell_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = drop.doublet)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "drop.doublet__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_B)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Inhib_B_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_H)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_H_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = drop.lowNTx_B)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "drop.lowNTx_B__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = drop.lowNTx_A)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "drop.lowNTx_A__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_G)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_G_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Mural)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Mural_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_F)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_F_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = OPC_COP)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "OPC_COP_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_C)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Inhib_C_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_A)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_A_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Inhib_D)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Inhib_D_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_E)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_E_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = Excit_D)) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Excit_D_cell_type__proportions") +
    facet_wrap(vars(sample_id)) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
