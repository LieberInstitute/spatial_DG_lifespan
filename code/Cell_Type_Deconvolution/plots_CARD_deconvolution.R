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
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_CARD_sestan.rds"))

#############################################################
# Scatterpie plots for deconvoluted cell types of Visium data
#############################################################

Br1412 <- spe[, spe$sample_id == "Br1412"]
Br1412y <- as.matrix(colData(Br1412)[, c(44:70)])

Br2706 <- spe[, spe$sample_id == "Br2706"]
Br2706y <- as.matrix(colData(Br2706)[, c(44:70)])

Br2720 <- spe[, spe$sample_id == "Br2720"]
Br2720y <- as.matrix(colData(Br2720)[, c(44:70)])

Br3942 <- spe[, spe$sample_id == "Br3942"]
Br3942y <- as.matrix(colData(Br3942)[, c(44:70)])

Br5242 <- spe[, spe$sample_id == "Br5242"]
Br5242y <- as.matrix(colData(Br5242)[, c(44:70)])

Br5699_new <- spe[, spe$sample_id == "Br5699_new"]
Br5699_newy <- as.matrix(colData(Br5699_new)[, c(44:70)])

Br6023 <- spe[, spe$sample_id == "Br6023"]
Br6023y <- as.matrix(colData(Br6023)[, c(44:70)])

Br6129_new <- spe[, spe$sample_id == "Br6129_new"]
Br6129_newy <- as.matrix(colData(Br6129_new)[, c(44:70)])

Br6299_new <- spe[, spe$sample_id == "Br6299_new"]
Br6299_newy <- as.matrix(colData(Br6299_new)[, c(44:70)])

Br6522 <- spe[, spe$sample_id == "Br6522"]
Br6522y <- as.matrix(colData(Br6522)[, c(44:70)])

Br8181 <- spe[, spe$sample_id == "Br8181"]
Br8181y <- as.matrix(colData(Br8181)[, c(44:70)])

Br8195 <- spe[, spe$sample_id == "Br8195"]
Br8195y <- as.matrix(colData(Br8195)[, c(44:70)])

Br8533 <- spe[, spe$sample_id == "Br8533"]
Br8533y <- as.matrix(colData(Br8533)[, c(44:70)])

Br8667 <- spe[, spe$sample_id == "Br8667"]
Br8667y <- as.matrix(colData(Br8667)[, c(44:70)])

Br8686 <- spe[, spe$sample_id == "Br8686"]
Br8686y <- as.matrix(colData(Br8686)[, c(44:70)])

Br8700 <- spe[, spe$sample_id == "Br8700"]
Br8700y <- as.matrix(colData(Br8700)[, c(44:70)])


cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3", "InN_NPY" = "green3",
    "VLMC" = "firebrick1", "aSMC" = "red1","Pericyte" = "red2", "Endoth" = "red", "vSMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "yellow", "Astro_2" = "yellow3")


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

#############################################################
# Plot cell-type proportions separately for each capture area
#############################################################

spe_cells <- as.data.frame(cbind(colData(spe)$sample_id, colData(spe)[, c(44:70)], spatialCoords(spe)))
spe_cells <- rename(spe_cells, "sample_id" = "colData.spe..sample_id")

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Cell_type_full_scale_proportions.pdf"))

ggplot(spe_cells) +
    geom_point(aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
        color = CA3_N), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "CA3_N_cell_type__proportions") +
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
        color = CA1_N), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "CA1_N_cell_type__proportions") +
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
        color = Mossy), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Mossy_cell_type__proportions") +
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
        color = CA2_N), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "CA2_N_cell_type__proportions") +
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
        color = GC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "GC_cell_type__proportions") +
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
        color = InN_LAMP5), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_LAMP5_cell_type__proportions") +
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
        color = InN_PV), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_PV_cell_type__proportions") +
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
        color = InN_LHX6), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_LHX6_cell_type__proportions") +
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
        color = InN_VIP), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_VIP_cell_type__proportions") +
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
        color = InN_SST), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_SST_cell_type__proportions") +
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
        color = InN_NR2F2), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_NR2F2_cell_type__proportions") +
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
        color = InN_MEIS2), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_MEIS2_cell_type__proportions") +
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
        color = InN_NPY), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "InN_NPY_cell_type__proportions") +
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
        color = Oligo), size = 0.095) +
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
        color = Microglia), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Microglia_cell_type__proportions") +
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
        color = OPC), size = 0.095) +
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
        color = Endoth), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Endoth_cell_type__proportions") +
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
        color = vSMC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "vSMC_Muscle_cell_type__proportions") +
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
        color = aSMC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "aSMC_Muscle_cell_type__proportions") +
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
        color = VLMC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "VLMC_cell_type__proportions") +
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
        color = Macro), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Macro_cell_type__proportions") +
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
        color = COP), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "COP_cell_type__proportions") +
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
        color = T_Cell), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "T_Cell_type__proportions") +
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
        color = Pericyte), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Pericyte_cell_type__proportions") +
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
        color = Myeloid), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Myeloid_cell_type__proportions") +
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
        color = Astro_1), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Astro_1_cell_type__proportions") +
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
        color = Astro_2), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion", limits = c(0, 1)) +
    labs(title = "Astro_2_cell_type__proportions") +
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
        color = CA3_N), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "CA3_N_cell_type__proportions") +
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
        color = CA1_N), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "CA1_N_cell_type__proportions") +
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
        color = Mossy), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Mossy_cell_type__proportions") +
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
        color = CA2_N), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "CA2_N_cell_type__proportions") +
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
        color = GC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "GC_cell_type__proportions") +
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
        color = InN_LAMP5), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_LAMP5_cell_type__proportions") +
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
        color = InN_PV), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_PV_cell_type__proportions") +
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
        color = InN_LHX6), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_LHX6_cell_type__proportions") +
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
        color = InN_VIP), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_VIP_cell_type__proportions") +
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
        color = InN_SST), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_SST_cell_type__proportions") +
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
        color = InN_NR2F2), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_NR2F2_cell_type__proportions") +
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
        color = InN_MEIS2), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_MEIS2_cell_type__proportions") +
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
        color = InN_NPY), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "InN_NPY_cell_type__proportions") +
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
        color = Oligo), size = 0.095) +
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
        color = Microglia), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion",) +
    labs(title = "Microglia_cell_type__proportions") +
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
        color = OPC), size = 0.095) +
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
        color = Endoth), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Endoth_cell_type__proportions") +
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
        color = vSMC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "vSMC_Muscle_cell_type__proportions") +
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
        color = aSMC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "aSMC_Muscle_cell_type__proportions") +
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
        color = VLMC), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "VLMC_cell_type__proportions") +
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
        color = Macro), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Macro_cell_type__proportions") +
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
        color = COP), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "COP_cell_type__proportions") +
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
        color = T_Cell), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "T_Cell_type__proportions") +
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
        color = Pericyte), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Pericyte_cell_type__proportions") +
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
        color = Myeloid), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Myeloid_cell_type__proportions") +
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
        color = Astro_1), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Astro_1_cell_type__proportions") +
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
        color = Astro_2), size = 0.095) +
    scale_color_viridis(option = "turbo", name = "proportion") +
    labs(title = "Astro_2_cell_type__proportions") +
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
