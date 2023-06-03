#########################################
# spatial_DG_lifespan project
# Radial distances of genes and celltypes
# Anthony Ramnauth, May 31 2023
#########################################

suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(Seurat)
    library(SeuratData)
    library(sessioninfo)
    library(semla)
    library(tibble)
    library(ggplot2)
    library(patchwork)
    library(scico)
    library(tidyr)
    library(dplyr)
})

# Load SPE
spe <-
    readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

# Add row and col columns in colData(spe), specific names needed for meta.data
spe$row <- spe$array_row
spe$col <- spe$array_col

# To make sure barcodes are unique use spe$key
colnames(spe) <- spe$key

DG_seu <- CreateSeuratObject(
    counts = as.matrix(counts(spe)),
    logcounts = as.matrix(logcounts(spe)),
    meta.data = data.frame(colData(spe)),
    feature.data = data.frame(rowData(spe)),
    img.data = data.frame(imgData(spe)),
    project = "DG_lifespan"
    )

seuList <- list()

samples <- unique(spe$sample_id)

for (i in seq_along(samples)){
seuList[[i]] = subset(x = DG_seu, subset = sample_id == samples[i])
}

names(seuList) <- samples

save(seuList, file = here("processed-data", "Seurat", "seuList_cell2loc.Rdata"))

# Format objects for calculating radial distances

for (s in seq_along(samples)){
u <- tibble(
    barcode = colnames(seuList[[s]]),
    x = as.numeric(seuList[[s]]$col),
    y = as.numeric(seuList[[s]]$row),
    sampleID = as.numeric(seuList[[s]]$age)
)

v <- seuList[[s]]$key[seuList[[s]]$bayesSpace_harmony_10 == "7"]

# Calculate radial distances
seuList[[s]]$radial_distances <- RadialDistance(u, v, convert_to_microns = TRUE)

}

# Transfer results to spe object
total_meta.data <- as.data.frame(rbind(seuList[[1]]@meta.data,
    seuList[[2]]@meta.data,
    seuList[[3]]@meta.data,
    seuList[[4]]@meta.data,
    seuList[[5]]@meta.data,
    seuList[[6]]@meta.data,
    seuList[[7]]@meta.data,
    seuList[[8]]@meta.data,
    seuList[[9]]@meta.data,
    seuList[[10]]@meta.data,
    seuList[[11]]@meta.data,
    seuList[[12]]@meta.data,
    seuList[[13]]@meta.data,
    seuList[[14]]@meta.data,
    seuList[[15]]@meta.data,
    seuList[[16]]@meta.data
        ))

length(setdiff(spe$key, total_meta.data$key))
setdiff(spe$key, total_meta.data$key)
spe$radial_distances <- total_meta.data$radial_distances

spe$r_dist_sqrt <- sign(spe$radial_distances)*sqrt(abs(spe$radial_distances))

vis_grid_gene(spe = spe,
    geneid = "radial_distances",
    spatial = TRUE,
    viridis - FALSE,
    pdf = here::here("plots", "spatstats", "GCL_radial_distance.pdf"),
    cont_colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"),
    image_id = "lowres",
    alpha = 0.5,
    minCount = -50,
    point_size = 2,
    sample_order = unique(spe$sample_id)
)

# Save spe object
saveRDS(spe, file = here("processed-data", "semla", "spe_radial_dist.rds"))

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

df <- as.data.frame(colData(spe)[,44:73])

# Plot the radial distances from GCL boarders of all mean cell abundances faceted by age_bin
pdf(file = here::here("plots", "spatstats", "Radial_dist_GCL_ML_celltypes.pdf"))

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Astro_1, color = "Astro_1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Astro_2, color = "Astro_2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_CA1_N, color = "CA1_N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_CA2_N, color = "CA2_N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_CA3_N, color = "CA3_N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_COP, color = "COP"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Endoth, color = "Endoth"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_GC, color = "GC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_LAMP5, color = "InN_LAMP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_LHX6, color = "InN_LHX6"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_MEIS2, color = "InN_MEIS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_NR2F2, color = "InN_NR2F2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_PV, color = "InN_PV"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_SST, color = "InN_SST"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_VIP, color = "InN_VIP"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Macro, color = "Macro"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Microglia, color = "Microglia"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Mossy, color = "Mossy"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Myeloid, color = "Myeloid"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_OPC, color = "OPC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Oligo, color = "Oligo"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Pericyte, color = "Pericyte"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_SMC, color = "SMC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_T_Cell, color = "T_Cell"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_VLMC, color = "VLMC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-500, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Astro_1, color = "Astro_1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Astro_2, color = "Astro_2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_CA1_N, color = "CA1_N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_CA2_N, color = "CA2_N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_CA3_N, color = "CA3_N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_COP, color = "COP"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Endoth, color = "Endoth"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_GC, color = "GC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_LAMP5, color = "InN_LAMP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_LHX6, color = "InN_LHX6"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_MEIS2, color = "InN_MEIS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_NR2F2, color = "InN_NR2F2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_PV, color = "InN_PV"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_SST, color = "InN_SST"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_InN_VIP, color = "InN_VIP"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Macro, color = "Macro"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Microglia, color = "Microglia"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Mossy, color = "Mossy"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Myeloid, color = "Myeloid"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_OPC, color = "OPC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Oligo, color = "Oligo"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_Pericyte, color = "Pericyte"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_SMC, color = "SMC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_T_Cell, color = "T_Cell"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

ggplot(df) +
    geom_smooth(aes(radial_distances, meanscell_abundance_w_sf_VLMC, color = "VLMC"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 4000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean cell abundance", x= "Radial distance from GCL boarder",
        title = "Radial distances of mean cell abundances from GCL")

dev.off()

################################################################################

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

# Create separate spe object for  BayesSpace spatial domain of interest
spe_ML_GCL <- spe[, which(spe$bayesSpace_harmony_10 == "2" | spe$bayesSpace_harmony_10 == "7")]
dim(spe_ML_GCL)

# Get list of gene-set from mouse data (Robert R. Stickels et al., 2020) for dendritically enriched gene sets
Stickels_2020 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Stickels_2020.csv"))

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Translate from one species to the other using the orthology
Dcluster_1 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.1,]
Dcluster_2 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.2,]
Dcluster_3 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.3,]
Dcluster_4 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster4,]


#############################################################################################
# Subset for dendritic cluster 3
spe_ML_d3 <- spe_ML_GCL[which((rownames(spe_ML_GCL)) %in% Dcluster_3$Column1)]

df_rad_ML_d3 <- as.data.frame(colData(spe_ML_d3)[,71:73])

# Make dataframe for all the logcounts and radial distances
for (g in seq_along(rownames(spe_ML_d3))) {
    Dcluster_3_df <- cbind(df_rad_ML_d3, t(as.data.frame(logcounts(spe_ML_d3))))
}

# Plot for all the cluster 3 genes faceted by age_bin
pdf(file = here::here("plots", "spatstats", "Radial_dist_GCL_ML_dendritic3.pdf"))

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, ANAPC13, color = "ANAPC13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CAMK2A, color = "CAMK2A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TULP4, color = "TULP4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NCS1, color = "NCS1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NSMF, color = "NSMF"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ZNF365, color = "ZNF365"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PSD, color = "PSD"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EIF3F, color = "EIF3F"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RDX, color = "RDX"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DDN, color = "DDN"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, KIF5A, color = "KIF5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, HOMER2, color = "HOMER2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL26, color = "RPL26"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PPP1R9B, color = "PPP1R9B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CSNK1G2, color = "CSNK1G2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF2, color = "EEF2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SLC25A23, color = "SLC25A23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SIRT2, color = "SIRT2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DACT3, color = "DACT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CNOT3, color = "CNOT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CENPB, color = "CENPB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NNAT, color = "NNAT"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PHACTR3, color = "PHACTR3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CABP7, color = "CABP7"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, WDR13, color = "WDR13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD2, color = "GLUD2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

# Without legend

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, ANAPC13, color = "ANAPC13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CAMK2A, color = "CAMK2A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TULP4, color = "TULP4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NCS1, color = "NCS1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NSMF, color = "NSMF"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ZNF365, color = "ZNF365"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PSD, color = "PSD"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EIF3F, color = "EIF3F"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RDX, color = "RDX"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DDN, color = "DDN"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, KIF5A, color = "KIF5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, HOMER2, color = "HOMER2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL26, color = "RPL26"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PPP1R9B, color = "PPP1R9B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CSNK1G2, color = "CSNK1G2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF2, color = "EEF2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SLC25A23, color = "SLC25A23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SIRT2, color = "SIRT2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DACT3, color = "DACT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CNOT3, color = "CNOT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CENPB, color = "CENPB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NNAT, color = "NNAT"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PHACTR3, color = "PHACTR3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CABP7, color = "CABP7"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, WDR13, color = "WDR13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD2, color = "GLUD2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    theme(legend.position="none") +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

#Plot each individual gene from dendritic cluster 3

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, ANAPC13, color = "ANAPC13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, CAMK2A, color = "CAMK2A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, TULP4, color = "TULP4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, NCS1, color = "NCS1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, NSMF, color = "NSMF"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, ZNF365, color = "ZNF365"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, PSD, color = "PSD"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, EIF3F, color = "EIF3F"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, RDX, color = "RDX"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, DDN, color = "DDN"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, KIF5A, color = "KIF5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, HOMER2, color = "HOMER2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, RPL26, color = "RPL26"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, PPP1R9B, color = "PPP1R9B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, CSNK1G2, color = "CSNK1G2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, EEF2, color = "EEF2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, SLC25A23, color = "SLC25A23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, SIRT2, color = "SIRT2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, DACT3, color = "DACT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, CNOT3, color = "CNOT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, CENPB, color = "CENPB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, NNAT, color = "NNAT"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, PHACTR3, color = "PHACTR3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, CABP7, color = "CABP7"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, WDR13, color = "WDR13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_3_df) +
    geom_smooth(aes(radial_distances, GLUD2, color = "GLUD2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 3 dendritically enriched gene set in GCL & ML")

dev.off()

###########################################################################################

# Subset for dendritic cluster 4
spe_ML_d4 <- spe_ML_GCL[which((rownames(spe_ML_GCL)) %in% Dcluster_4$Column1)]

df_rad_ML_d4 <- as.data.frame(colData(spe_ML_d4)[,71:73])

# Make dataframe for all the logcounts and radial distances
for (g in seq_along(rownames(spe_ML_d4))) {
    Dcluster_4_df <- cbind(df_rad_ML_d4, t(as.data.frame(logcounts(spe_ML_d4))))
}

# Plot for all the cluster 3 genes faceted by age_bin
pdf(file = here::here("plots", "spatstats", "Radial_dist_GCL_ML_dendritic4.pdf"))

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPL22, color = "RPL22"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CAMK2N1, color = "CAMK2N1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TRNP1, color = "TRNP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, HPCAL4, color = "HPCAL4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GAS5, color = "GAS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, AKT3, color = "AKT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, UGP2, color = "UGP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP4K4, color = "MAP4K4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP2, color = "MAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ACSL3, color = "ACSL3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CHST2, color = "CHST2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PPP2R2C, color = "PPP2R2C"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3A, color = "RPS3A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SOWAHA, color = "SOWAHA"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NRSN1, color = "NRSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PURB, color = "PURB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DENND5A, color = "DENND5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAPK8IP1, color = "MAPK8IP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FTH1, color = "FTH1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX8A, color = "COX8A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CNIH2, color = "CNIH2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PITPNM1, color = "PITPNM1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CACNB3, color = "CACNB3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, AGAP2, color = "AGAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, WSB2, color = "WSB2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX6A1, color = "COX6A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RNF10, color = "RNF10"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, JPH4, color = "JPH4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CEP170B, color = "CEP170B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ZNF106, color = "ZNF106"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL4, color = "RPL4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FBXL16, color = "FBXL16"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ST3GAL2, color = "ST3GAL2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ZBTB4, color = "ZBTB4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GIT1, color = "GIT1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL23, color = "RPL23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ABHD17A, color = "ABHD17A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ABHD8, color = "ABHD8"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SNPH, color = "SNPH"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DESI1, color = "DESI1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PCSK1N, color = "PCSK1N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, IDS, color = "IDS"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

# Without legend

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPL22, color = "RPL22"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CAMK2N1, color = "CAMK2N1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TRNP1, color = "TRNP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, HPCAL4, color = "HPCAL4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GAS5, color = "GAS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, AKT3, color = "AKT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, UGP2, color = "UGP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP4K4, color = "MAP4K4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP2, color = "MAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ACSL3, color = "ACSL3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CHST2, color = "CHST2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PPP2R2C, color = "PPP2R2C"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3A, color = "RPS3A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SOWAHA, color = "SOWAHA"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NRSN1, color = "NRSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PURB, color = "PURB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DENND5A, color = "DENND5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAPK8IP1, color = "MAPK8IP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FTH1, color = "FTH1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX8A, color = "COX8A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CNIH2, color = "CNIH2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PITPNM1, color = "PITPNM1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CACNB3, color = "CACNB3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, AGAP2, color = "AGAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, WSB2, color = "WSB2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX6A1, color = "COX6A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RNF10, color = "RNF10"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, JPH4, color = "JPH4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CEP170B, color = "CEP170B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ZNF106, color = "ZNF106"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL4, color = "RPL4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FBXL16, color = "FBXL16"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ST3GAL2, color = "ST3GAL2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ZBTB4, color = "ZBTB4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GIT1, color = "GIT1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL23, color = "RPL23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ABHD17A, color = "ABHD17A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ABHD8, color = "ABHD8"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SNPH, color = "SNPH"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DESI1, color = "DESI1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PCSK1N, color = "PCSK1N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, IDS, color = "IDS"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    theme(legend.position="none") +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

#Plot each individual gene from dendritic cluster 4

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPL22, color = "RPL22"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, CAMK2N1, color = "CAMK2N1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, TRNP1, color = "TRNP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, HPCAL4, color = "HPCAL4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, GAS5, color = "GAS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, AKT3, color = "AKT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, UGP2, color = "UGP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, MAP4K4, color = "MAP4K4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, MAP2, color = "MAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ACSL3, color = "ACSL3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, CHST2, color = "CHST2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, PPP2R2C, color = "PPP2R2C"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS3A, color = "RPS3A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, SOWAHA, color = "SOWAHA"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, NRSN1, color = "NRSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, PURB, color = "PURB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, DENND5A, color = "DENND5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, MAPK8IP1, color = "MAPK8IP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, FTH1, color = "FTH1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, COX8A, color = "COX8A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, CNIH2, color = "CNIH2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, PITPNM1, color = "PITPNM1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, CACNB3, color = "CACNB3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, AGAP2, color = "AGAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, WSB2, color = "WSB2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, COX6A1, color = "COX6A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RNF10, color = "RNF10"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, JPH4, color = "JPH4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, CEP170B, color = "CEP170B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ZNF106, color = "ZNF106"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPL4, color = "RPL4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, FBXL16, color = "FBXL16"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ST3GAL2, color = "ST3GAL2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ZBTB4, color = "ZBTB4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, GIT1, color = "GIT1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPL23, color = "RPL23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ABHD17A, color = "ABHD17A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ABHD8, color = "ABHD8"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, SNPH, color = "SNPH"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, DESI1, color = "DESI1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, PCSK1N, color = "PCSK1N"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

ggplot(Dcluster_4_df) +
    geom_smooth(aes(radial_distances, IDS, color = "IDS"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Cluster 4 dendritically enriched gene set in GCL & ML")

dev.off()

##########################################################################################

# Isolate the dendritic markers that are very enriched from the above plots in ML

combined_gene_set <- c(
    "CAMK2A", "TUBB2B", "EEF1A1", "TULP4", "LHX2", "NCS1", "NSMF", "GLUD1", "PSD",
    "RDX", "DDN", "ARHGAP5", "HOMER2", "RPL13", "RPL26", "PPP1R9B", "RPS15", "EEF2",
    "SLC25A23", "SIRT2", "DACT3", "RPL13A", "CENPB", "CABP7", "RPL22", "MAP4K4",
    "MAP2", "CCNI", "RPS3A", "RPS14", "RPS24", "DENND5A", "MAPK8IP1", "FTH1", "EEF1G",
    "COX8A", "PITPNM1", "RPS3", "AGAP2", "RPLP0", "COX6A1", "CEP170B", "RPLP1", "FBXL16",
    "RPS2", "COX4I1", "RPL23", "ABHD17A", "RPS11", "RPS9", "RPS5", "ITSN1", "GPM6B",
    "IDS"
)

spe_ML_dendcombined <- spe_ML_GCL[which((rownames(spe_ML_GCL)) %in% combined_gene_set)]
df_rad_ML_dendcombined <- as.data.frame(colData(spe_ML_dendcombined)[,71:73])

# Make dataframe for all the logcounts and radial distances
for (g in seq_along(rownames(spe_ML_dendcombined))) {
    combined_dend_df <- cbind(df_rad_ML_dendcombined, t(as.data.frame(logcounts(spe_ML_dendcombined))))
}

# Plot for all the cluster 3 genes faceted by age_bin
pdf(file = here::here("plots", "spatstats", "Radial_dist_GCL_ML_dendritic_select.pdf"))

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL22, color = "RPL22"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP4K4, color = "MAP4K4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP2, color = "MAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3A, color = "RPS3A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CAMK2A, color = "CAMK2A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TULP4, color = "TULP4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NCS1, color = "NCS1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NSMF, color = "NSMF"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PSD, color = "PSD"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DENND5A, color = "DENND5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAPK8IP1, color = "MAPK8IP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FTH1, color = "FTH1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX8A, color = "COX8A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PITPNM1, color = "PITPNM1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RDX, color = "RDX"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DDN, color = "DDN"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, AGAP2, color = "AGAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX6A1, color = "COX6A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CEP170B, color = "CEP170B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, HOMER2, color = "HOMER2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FBXL16, color = "FBXL16"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL26, color = "RPL26"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL23, color = "RPL23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PPP1R9B, color = "PPP1R9B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ABHD17A, color = "ABHD17A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF2, color = "EEF2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SLC25A23, color = "SLC25A23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SIRT2, color = "SIRT2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DACT3, color = "DACT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CENPB, color = "CENPB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CABP7, color = "CABP7"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, IDS, color = "IDS"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

# Without legend

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL22, color = "RPL22"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP4K4, color = "MAP4K4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAP2, color = "MAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3A, color = "RPS3A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CAMK2A, color = "CAMK2A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, TULP4, color = "TULP4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NCS1, color = "NCS1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, NSMF, color = "NSMF"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PSD, color = "PSD"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DENND5A, color = "DENND5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, MAPK8IP1, color = "MAPK8IP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FTH1, color = "FTH1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX8A, color = "COX8A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PITPNM1, color = "PITPNM1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RDX, color = "RDX"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DDN, color = "DDN"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, AGAP2, color = "AGAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX6A1, color = "COX6A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CEP170B, color = "CEP170B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, HOMER2, color = "HOMER2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, FBXL16, color = "FBXL16"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL26, color = "RPL26"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL23, color = "RPL23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, PPP1R9B, color = "PPP1R9B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ABHD17A, color = "ABHD17A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF2, color = "EEF2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SLC25A23, color = "SLC25A23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, SIRT2, color = "SIRT2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, DACT3, color = "DACT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CENPB, color = "CENPB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CABP7, color = "CABP7"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, IDS, color = "IDS"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    theme(legend.position="none") +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

#Plot each individual gene from dendritic cluster 4

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL22, color = "RPL22"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, MAP4K4, color = "MAP4K4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, MAP2, color = "MAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS3A, color = "RPS3A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, CAMK2A, color = "CAMK2A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, TULP4, color = "TULP4"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, NCS1, color = "NCS1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, NSMF, color = "NSMF"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, PSD, color = "PSD"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, DENND5A, color = "DENND5A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, MAPK8IP1, color = "MAPK8IP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, FTH1, color = "FTH1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, COX8A, color = "COX8A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, PITPNM1, color = "PITPNM1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RDX, color = "RDX"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, DDN, color = "DDN"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, AGAP2, color = "AGAP2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, COX6A1, color = "COX6A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, CEP170B, color = "CEP170B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, HOMER2, color = "HOMER2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, FBXL16, color = "FBXL16"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL26, color = "RPL26"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL23, color = "RPL23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, PPP1R9B, color = "PPP1R9B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, ABHD17A, color = "ABHD17A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, EEF2, color = "EEF2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, SLC25A23, color = "SLC25A23"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, SIRT2, color = "SIRT2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, DACT3, color = "DACT3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL & ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, CENPB, color = "CENPB"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, CABP7, color = "CABP7"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

ggplot(combined_dend_df) +
    geom_smooth(aes(radial_distances, IDS, color = "IDS"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically enriched gene set in GCL to ML")

dev.off()

##################################################################################################################

# Look at elderly depleted genes

# Isolate the dendritic markers that are very enriched from the above plots in ML

depleted_gene_set <- c(
    "TUBB2B", "EEF1A1", "LHX2", "GLUD1", "ARHGAP5", "RPL13", "RPS15",
    "RPL13A", "CCNI", "RPS14", "RPS24", "EEF1G", "RPS3", "RPLP0", "RPLP1",
    "RPS2", "COX4I1", "RPS11", "RPS9", "RPS5", "ITSN1", "GPM6B"
)

spe_ML_depleted <- spe_ML_GCL[which((rownames(spe_ML_GCL)) %in% depleted_gene_set)]
df_rad_ML_depleted <- as.data.frame(colData(spe_ML_depleted)[,71:73])

# Make dataframe for all the logcounts and radial distances
for (g in seq_along(rownames(spe_ML_depleted))) {
    depleted_dend_df <- cbind(df_rad_ML_depleted, t(as.data.frame(logcounts(spe_ML_depleted))))
}

# Plot for all the cluster 3 genes faceted by age_bin
pdf(file = here::here("plots", "spatstats", "Radial_dist_GCL_ML_dendritic_depleted.pdf"))

ggplot(depleted_dend_df) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically age-depleting gene set in GCL to ML")

ggplot(depleted_dend_df) +
    geom_smooth(aes(radial_distances, TUBB2B, color = "TUBB2B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1A1, color = "EEF1A1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, LHX2, color = "LHX2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GLUD1, color = "GLUD1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ARHGAP5, color = "ARHGAP5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13, color = "RPL13"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS15, color = "RPS15"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPL13A, color = "RPL13A"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, CCNI, color = "CCNI"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS14, color = "RPS14"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS24, color = "RPS24"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, EEF1G, color = "EEF1G"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS3, color = "RPS3"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP0, color = "RPLP0"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPLP1, color = "RPLP1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS2, color = "RPS2"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, COX4I1, color = "COX4I1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS11, color = "RPS11"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS9, color = "RPS9"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, RPS5, color = "RPS5"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, ITSN1, color = "ITSN1"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(aes(radial_distances, GPM6B, color = "GPM6B"),
        method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 0), linetype = "dashed") +
    xlim(-750, 1000) +
    facet_wrap(vars(age_bin)) +
    theme_bw() +
    theme_classic() +
    theme(legend.position="none") +
    labs(y = "Mean normalize logcount", x= "Radial distance from GCL boarder",
        title = "Select conserved dendritically age-depleting gene set in GCL to ML")

dev.off()

