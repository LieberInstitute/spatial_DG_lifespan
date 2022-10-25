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
    library(ggcorrplot)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

#############################################################
# Scatterpie plots for deconvoluted cell types of Visium data
#############################################################

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

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Scatterpie_plots.pdf"), width = 12, height = 8)
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
        color = Oligo), size = 0.09) +
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
        color = Micro), size = 0.09) +
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
        color = OPC), size = 0.09) +
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
        color = Inhib_A), size = 0.09) +
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
        color = Astro_A), size = 0.09) +
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
        color = Excit_B), size = 0.09) +
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
        color = Astro_B), size = 0.09) +
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
        color = Excit_C), size = 0.09) +
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
        color = Tcell), size = 0.09) +
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
        color = drop.doublet), size = 0.09) +
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
        color = Inhib_B), size = 0.09) +
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
        color = Excit_H), size = 0.09) +
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
        color = drop.lowNTx_B), size = 0.09) +
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
        color = drop.lowNTx_A), size = 0.09) +
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
        color = Excit_G), size = 0.09) +
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
        color = Mural), size = 0.09) +
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
        color = Excit_F), size = 0.09) +
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
        color = OPC_COP), size = 0.09) +
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
        color = Inhib_C), size = 0.09) +
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
        color = Excit_A), size = 0.09) +
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
        color = Inhib_D), size = 0.09) +
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
        color = Excit_E), size = 0.09) +
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
        color = Excit_D), size = 0.09) +
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
        color = Oligo), size = 0.09) +
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
        color = Micro), size = 0.09) +
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
        color = OPC), size = 0.09) +
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
        color = Inhib_A), size = 0.09) +
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
        color = Astro_A), size = 0.09) +
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
        color = Excit_B), size = 0.09) +
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
        color = Astro_B), size = 0.09) +
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
        color = Excit_C), size = 0.09) +
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
        color = Tcell), size = 0.09) +
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
        color = drop.doublet), size = 0.09) +
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
        color = Inhib_B), size = 0.09) +
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
        color = Excit_H), size = 0.09) +
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
        color = drop.lowNTx_B), size = 0.09) +
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
        color = drop.lowNTx_A), size = 0.09) +
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
        color = Excit_G), size = 0.09) +
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
        color = Mural), size = 0.09) +
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
        color = Excit_F), size = 0.09) +
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
        color = OPC_COP), size = 0.09) +
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
        color = Inhib_C), size = 0.09) +
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
        color = Excit_A), size = 0.09) +
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
        color = Inhib_D), size = 0.09) +
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
        color = Excit_E), size = 0.09) +
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
        color = Excit_D), size = 0.09) +
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

############################################################################
# Plot correlation matrix for cell-type proportions across spatial locations
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

all_spe_cells <- spe_cells
all_spe_cells$sample_id <- NULL
all_spe_cells$pxl_col_in_fullres <- NULL
all_spe_cells$pxl_row_in_fullres <- NULL

#reorder correlation by column index
all_spe_cells <- all_spe_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]

cor_all_spe_cells = cor(as.matrix(all_spe_cells))

# Age groups for all Visium spots
infant_spe <- spe[, spe$age_bin == "Infant"]
infant_cells <- as.data.frame(colData(infant_spe)[, c(44:66)])
infant_cells <- infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_infant_cells = cor(as.matrix(infant_cells))

teen_spe <- spe[, spe$age_bin == "Teen"]
teen_cells <- as.data.frame(colData(teen_spe)[, c(44:66)])
teen_cells <- teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_teen_cells = cor(as.matrix(teen_cells))

adult_spe <- spe[, spe$age_bin == "Adult"]
adult_cells <- as.data.frame(colData(adult_spe)[, c(44:66)])
adult_cells <- adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_adult_cells = cor(as.matrix(adult_cells))

elderly_spe <- spe[, spe$age_bin == "Elderly"]
elderly_cells <- as.data.frame(colData(elderly_spe)[, c(44:66)])
elderly_cells <- elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_elderly_cells = cor(as.matrix(elderly_cells))


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_correlation_matrix_all_spots.pdf"),
    width = 12, height = 8)
ggcorrplot(cor_all_spe_cells, hc.order = F,
    outline.color = "white",
    tl.srt = 60,
    tl.cex = 18,
    lab_size = 7,
    colors = c("#91a28c","white","#8f2c37"))+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 16),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm')) +
    coord_fixed()+
    ggtitle("Correlations in cell-type proportion across all Visium spatial locations") +
    theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_infant_cells, hc.order = F,
    outline.color = "white",
    tl.srt = 60,
    tl.cex = 18,
    lab_size = 7,
    colors = c("#91a28c","white","#8f2c37"))+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 16),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm')) +
    coord_fixed()+
    ggtitle("Correlations in cell-type proportion across Infant Visium spatial locations") +
    theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_teen_cells, hc.order = F,
    outline.color = "white",
    tl.srt = 60,
    tl.cex = 18,
    lab_size = 7,
    colors = c("#91a28c","white","#8f2c37"))+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 16),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm')) +
    coord_fixed()+
    ggtitle("Correlations in cell-type proportion across Teen Visium spatial locations") +
    theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_adult_cells, hc.order = F,
    outline.color = "white",
    tl.srt = 60,
    tl.cex = 18,
    lab_size = 7,
    colors = c("#91a28c","white","#8f2c37"))+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 16),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm')) +
    coord_fixed()+
    ggtitle("Correlations in cell-type proportion across Adult Visium spatial locations") +
    theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_elderly_cells, hc.order = F,
    outline.color = "white",
    tl.srt = 60,
    tl.cex = 18,
    lab_size = 7,
    colors = c("#91a28c","white","#8f2c37"))+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 16),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm')) +
    coord_fixed()+
    ggtitle("Correlations in cell-type proportion across Elderly Visium spatial locations") +
    theme(plot.title = element_text(size=14,face="bold"))
dev.off()

##############################################
# subset spe data based on BayesSpace clusters
##############################################

###########
# Cluster 4
###########

spe_cluster4 <- spe[, spe$bayesSpace_harmony_8 == 4]
cluster4_cells <- as.data.frame(colData(spe_cluster4)[, c(44:66)])
cluster4_cells <- cluster4_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_cells = cor(as.matrix(cluster4_cells))

# Age groups for cluster 4
infant_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Infant"]
cluster4_infant_cells <- as.data.frame(colData(infant_cluster4)[, c(44:66)])
cluster4_infant_cells <- cluster4_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_infant_cells = cor(as.matrix(cluster4_infant_cells))

teen_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Teen"]
cluster4_teen_cells <- as.data.frame(colData(teen_cluster4)[, c(44:66)])
cluster4_teen_cells <- cluster4_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_teen_cells = cor(as.matrix(cluster4_teen_cells))

adult_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Adult"]
cluster4_adult_cells <- as.data.frame(colData(adult_cluster4)[, c(44:66)])
cluster4_adult_cells <- cluster4_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_adult_cells = cor(as.matrix(cluster4_adult_cells))

elderly_cluster4 <- spe_cluster4[, spe_cluster4$age_bin == "Elderly"]
cluster4_elderly_cells <- as.data.frame(colData(elderly_cluster4)[, c(44:66)])
cluster4_elderly_cells <- cluster4_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster4_elderly_cells = cor(as.matrix(cluster4_elderly_cells))


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_correlation_matrix_cluster4_spots.pdf"),
    width = 12, height = 8)
ggcorrplot(cor_cluster4_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across All BayesSpace cluster_4 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster4_infant_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Infant BayesSpace cluster_4 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster4_teen_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Teen BayesSpace cluster_4 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster4_adult_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Adult BayesSpace cluster_4 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster4_elderly_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Elderly BayesSpace cluster_4 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))
dev.off()

############
# Cluster 1
###########

spe_cluster1 <- spe[, spe$bayesSpace_harmony_8 == 1]
cluster1_cells <- as.data.frame(colData(spe_cluster1)[, c(44:66)])
cluster1_cells <- cluster1_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_cells = cor(as.matrix(cluster1_cells))

# Age groups for cluster 1
infant_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Infant"]
cluster1_infant_cells <- as.data.frame(colData(infant_cluster1)[, c(44:66)])
cluster1_infant_cells <- cluster1_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_infant_cells = cor(as.matrix(cluster1_infant_cells))

teen_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Teen"]
cluster1_teen_cells <- as.data.frame(colData(teen_cluster1)[, c(44:66)])
cluster1_teen_cells <- cluster1_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_teen_cells = cor(as.matrix(cluster1_teen_cells))

adult_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Adult"]
cluster1_adult_cells <- as.data.frame(colData(adult_cluster1)[, c(44:66)])
cluster1_adult_cells <- cluster1_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_adult_cells = cor(as.matrix(cluster1_adult_cells))

elderly_cluster1 <- spe_cluster1[, spe_cluster1$age_bin == "Elderly"]
cluster1_elderly_cells <- as.data.frame(colData(elderly_cluster1)[, c(44:66)])
cluster1_elderly_cells <- cluster1_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster1_elderly_cells = cor(as.matrix(cluster1_elderly_cells))


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_correlation_matrix_cluster1_spots.pdf"),
    width = 12, height = 8)
ggcorrplot(cor_cluster1_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across All BayesSpace cluster_1 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster1_infant_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Infant BayesSpace cluster_1 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster1_teen_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Teen BayesSpace cluster_1 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster1_adult_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Adult BayesSpace cluster_1 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster1_elderly_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Elderly BayesSpace cluster_1 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))
dev.off()

###########
# Cluster 2
###########

spe_cluster2 <- spe[, spe$bayesSpace_harmony_8 == 2]
cluster2_cells <- as.data.frame(colData(spe_cluster2)[, c(44:66)])
cluster2_cells <- cluster2_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_cells = cor(as.matrix(cluster2_cells))

# Age groups for cluster 2
infant_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Infant"]
cluster2_infant_cells <- as.data.frame(colData(infant_cluster2)[, c(44:66)])
cluster2_infant_cells <- cluster2_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_infant_cells = cor(as.matrix(cluster2_infant_cells))

teen_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Teen"]
cluster2_teen_cells <- as.data.frame(colData(teen_cluster2)[, c(44:66)])
cluster2_teen_cells <- cluster2_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_teen_cells = cor(as.matrix(cluster2_teen_cells))

adult_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Adult"]
cluster2_adult_cells <- as.data.frame(colData(adult_cluster2)[, c(44:66)])
cluster2_adult_cells <- cluster2_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_adult_cells = cor(as.matrix(cluster2_adult_cells))

elderly_cluster2 <- spe_cluster2[, spe_cluster2$age_bin == "Elderly"]
cluster2_elderly_cells <- as.data.frame(colData(elderly_cluster2)[, c(44:66)])
cluster2_elderly_cells <- cluster2_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster2_elderly_cells = cor(as.matrix(cluster2_elderly_cells))


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_correlation_matrix_cluster2_spots.pdf"),
    width = 12, height = 8)
ggcorrplot(cor_cluster2_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across All BayesSpace cluster_2 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster2_infant_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Infant BayesSpace cluster_2 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster2_teen_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Teen BayesSpace cluster_2 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster2_adult_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Adult BayesSpace cluster_2 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster2_elderly_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Elderly BayesSpace cluster_2 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))
dev.off()

###########
# Cluster 8
###########

spe_cluster8 <- spe[, spe$bayesSpace_harmony_8 == 8]
cluster8_cells <- as.data.frame(colData(spe_cluster8)[, c(44:66)])
cluster8_cells <- cluster8_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_cells = cor(as.matrix(cluster8_cells))

# Age groups for cluster 8
infant_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Infant"]
cluster8_infant_cells <- as.data.frame(colData(infant_cluster8)[, c(44:66)])
cluster8_infant_cells <- cluster8_infant_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_infant_cells = cor(as.matrix(cluster8_infant_cells))

teen_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Teen"]
cluster8_teen_cells <- as.data.frame(colData(teen_cluster8)[, c(44:66)])
cluster8_teen_cells <- cluster8_teen_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_teen_cells = cor(as.matrix(cluster8_teen_cells))

adult_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Adult"]
cluster8_adult_cells <- as.data.frame(colData(adult_cluster8)[, c(44:66)])
cluster8_adult_cells <- cluster8_adult_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_adult_cells = cor(as.matrix(cluster8_adult_cells))

elderly_cluster8 <- spe_cluster8[, spe_cluster8$age_bin == "Elderly"]
cluster8_elderly_cells <- as.data.frame(colData(elderly_cluster8)[, c(44:66)])
cluster8_elderly_cells <- cluster8_elderly_cells[, c("Oligo", "OPC", "OPC_COP", "Astro_A", "Astro_B", "Micro",
    "Tcell", "Mural", "Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F",
    "Excit_G", "Excit_H", "Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "drop.lowNTx_A",
    "drop.lowNTx_B", "drop.doublet")]
cor_cluster8_elderly_cells = cor(as.matrix(cluster8_elderly_cells))


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell_type_correlation_matrix_cluster8_spots.pdf"),
    width = 12, height = 8)
ggcorrplot(cor_cluster8_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across All BayesSpace cluster_8 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster8_infant_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Infant BayesSpace cluster_8 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster8_teen_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Teen BayesSpace cluster_8 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster8_adult_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Adult BayesSpace cluster_8 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))

ggcorrplot(cor_cluster8_elderly_cells, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = c("#91a28c","white","#8f2c37"))+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Correlations in cell-type proportion across Elderly BayesSpace cluster_8 spatial locations") +
 theme(plot.title = element_text(size=14,face="bold"))
dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
