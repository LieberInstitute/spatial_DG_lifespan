###################################################################
# spatial_DG_lifespan project
# Plotting cell2location results
# Anthony Ramnauth, May 01 2023
###################################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
library(SpatialExperiment)
library(here)
library(ggplot2)
library(viridis)
})

# directory to save plots
dir_plots <- here("plots", "Cell_Type_Deconvolution", "cell2location")


# ---------------
# load SPE object
# ---------------

# load SPE object containing cell2location results
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

# check
head(colData(spe), 1)

# samples
sample_ids <- unique(colData(spe)$sample_id)
sample_ids

# names of columns containing deconvolved cell types
cols <- c(
    "meanscell_abundance_w_sf_Astro_1",
    "meanscell_abundance_w_sf_Astro_2",
    "meanscell_abundance_w_sf_CA1_N",
    "meanscell_abundance_w_sf_CA2_N",
    "meanscell_abundance_w_sf_CA3_N",
    "meanscell_abundance_w_sf_COP",
    "meanscell_abundance_w_sf_Endoth",
    "meanscell_abundance_w_sf_GC",
    "meanscell_abundance_w_sf_InN_LAMP5",
    "meanscell_abundance_w_sf_InN_LHX6",
    "meanscell_abundance_w_sf_InN_MEIS2",
    "meanscell_abundance_w_sf_InN_NR2F2",
    "meanscell_abundance_w_sf_InN_PV",
    "meanscell_abundance_w_sf_InN_SST",
    "meanscell_abundance_w_sf_InN_VIP",
    "meanscell_abundance_w_sf_Macro",
    "meanscell_abundance_w_sf_Microglia",
    "meanscell_abundance_w_sf_Mossy",
    "meanscell_abundance_w_sf_Myeloid",
    "meanscell_abundance_w_sf_OPC",
    "meanscell_abundance_w_sf_Oligo",
    "meanscell_abundance_w_sf_Pericyte",
    "meanscell_abundance_w_sf_SMC",
    "meanscell_abundance_w_sf_T_Cell",
    "meanscell_abundance_w_sf_VLMC"
)


# --------------
# generate plots
# --------------

# plot each cell type in each Visium sample

for (s in seq_along(sample_ids)) {

  # select sample
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]

  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))

  for (q in seq_along(cols)) {
    p <- ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres",
                               color = cols[q])) +
      geom_point(size = 0.5) +
      coord_fixed() +
      scale_y_reverse() +
      scale_color_viridis(option = "magma", name = "abundance") +
      labs(title = gsub("^.*sf_", "", cols[q]),
           subtitle = sample_ids[s]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "gray80"),
            panel.grid = element_line(color = "gray80"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())

    if (!dir.exists(here(dir_plots, sample_ids[s]))) {
      dir.create(here(dir_plots, sample_ids[s]), recursive = TRUE)
    }
    fn <- here(dir_plots, sample_ids[s],
               paste0(sample_ids[s], "_", gsub("^.*sf_", "", cols[q])))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 4, height = 3)
    ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 3)
  }
}


# plot each cell type across samples

min_all <- min(as.matrix(colData(spe)[, cols]))
max_all <- max(as.matrix(colData(spe)[, cols]))

for (q in seq_along(cols)) {

  # select cell type
  df <- as.data.frame(cbind(
    colData(spe)[, c("sample_id", cols[q]), drop = FALSE],
    spatialCoords(spe)
  ))

  p <- ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres",
                             color = cols[q])) +
    facet_wrap(~ sample_id, nrow = 2, scales = "free") +
    geom_point(size = 0.25) +
    scale_y_reverse() +
    scale_color_viridis(option = "magma", name = "abundance", trans = "sqrt",
                        limits = c(0, max_all)) +
    labs(title = gsub("^.*sf_", "", cols[q])) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "gray80"),
          panel.grid = element_line(color = "gray80"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

  fn <- here(dir_plots, paste0("cell2location_merged_", gsub("^.*sf_", "", cols[q])))
  ggsave(paste0(fn, ".pdf"), plot = p, width = 7, height = 4.75)
  ggsave(paste0(fn, ".png"), plot = p, width = 7, height = 4.75)
}
