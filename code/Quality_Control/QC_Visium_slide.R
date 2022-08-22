##################################
# spatial_DG_lifespan project
# QC individual Visium Slide
# Anthony Ramnauth, August 22 2022
##################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(scater)
    library(ggplot2)
    library(ggnewscale)
    library(ggspavis)
    library(scuttle)
    library(scran)
    library(dplyr)
    library(sessioninfo)
})

# load saved SPE object
spe <- readRDS(here::here("processed-data", "02_build_spe", "spe.rds"))

# Add variable of slide to colData(spe)
slide_df <- data.frame(spe$key, spe$sample_id)
slide_df <- slide_df %>%
    mutate(slide = case_when(
        grepl("Br1412", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br2706", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br3942", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br5242", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br5699", slide_df$spe.sample_id) ~ "V11D01387",
        grepl("Br6023", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br6129", slide_df$spe.sample_id) ~ "V11D01387",
        grepl("Br6299", slide_df$spe.sample_id) ~ "V11D01387",
        grepl("Br8195", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br8667", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br8686", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br8851", slide_df$spe.sample_id) ~ "V11D01387",
    ))

colData(spe)$slide <- factor(slide_df$slide, levels = c("V10B01088", "V11L05334", "V11D01387"))

## subset spe data based on slide
third_spe <- spe[, spe$slide %in% c("V11D01387")]

vis_grid_clus(
    spe = third_spe,
    clustervar = "in_tissue",
    pdf = here::here("plots", "QC_plots", "V11D01387_in_tissue_grid.pdf"),
    sort_clust = FALSE,
    colors = c("TRUE" = "grey90", "FALSE" = "orange"),
    point_size = 2
)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(third_spe)$gene_name)
table(is_mito)
rowData(third_spe)$gene_name[is_mito]

# calculate per-spot QC metrics and store in colData
third_spe <- addPerCellQC(third_spe, subsets = list(mito = is_mito))
head(colData(third_spe))

# Create QC Metrics
qcstats <- perCellQCMetrics(third_spe, subsets = list(
    Mito = which(seqnames(third_spe) == "chrM")
))

qcfilter <- DataFrame(
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = third_spe$sample_id),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = third_spe$sample_id),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher", batch = third_spe$sample_id)
)

colSums(as.matrix(qcfilter))

qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent

# Add QC metrics to spe object
third_spe$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
third_spe$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
third_spe$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
third_spe$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

# plot histograms of QC metrics
pdf(file = here::here("plots", "QC_plots", "V11D01387_QC_histograms_allSpots.pdf"), width = 8, height = 2.5)
par(mfrow = c(1, 3))
hist(colData(third_spe)$sum_umi, xlab = "sum", main = "UMIs per spot all donors")
hist(colData(third_spe)$sum_gene, xlab = "detected", main = "Genes per spot all donors")
hist(colData(third_spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito all donors")
par(mfrow = c(1, 1))
dev.off()

# QC plot of tissue spots discarded
for (i in colnames(qcfilter)) {
    vis_grid_clus(
        spe = third_spe,
        clustervar = paste0("scran_", i),
        pdf = here::here("plots", "QC_plots", paste0("V11D01387_scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        point_size = 2
    )
}


## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
