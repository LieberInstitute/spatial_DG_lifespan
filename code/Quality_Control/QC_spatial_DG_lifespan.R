#################################
# spatial_DG_lifespan project
# Script for quality control (QC)
# Anthony Ramnauth, Apr 12 2022
#################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

library(SpatialExperiment)
library(here)
library(scater)
library(ggplot2)
library(ggnewscale)
library(ggspavis)
library(scuttle)
library(scran)

# load saved SPE object
spe <- readRDS(here::here("processed-data", "02_build_spe", "spe.rds"))

# Create directory for QC plots
dir_rdata <- here::here("plots", "QC_plots")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

dim(spe)

# check SPE object contains only spots over tissue
table(colData(spe)$in_tissue)
all(colData(spe)$in_tissue)
dim(spe)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# Add Mito QC metrics using scater package
spe <- addPerCellQCMetrics(spe, subsets = list(mito = is_mito))

# Plot spatial coordinates of within tissue spots
pdf(file = here::here("plots", "QC_plots", "spe_Spatial_Coordinates.pdf"))
plotSpots(spe = spe[, colData(spe)$sample_id == "Br8686"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue") +
    ggtitle("Br8686")
plotSpots(spe = spe[, colData(spe)$sample_id == "Br2706"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue") +
    ggtitle("Br2706")
plotSpots(spe = spe[, colData(spe)$sample_id == "Br3942"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue") +
    ggtitle("Br3942")
plotSpots(spe = spe[, colData(spe)$sample_id == "Br6023"],
    x_coord = "pxl_col_in_fullres",
    y_coord = "pxl_row_in_fullres",
    size = 1,
    in_tissue = "in_tissue") +
    ggtitle("Br6023")
dev.off()

# plot histograms of QC metrics
pdf(file = here::here("plots", "QC_plots", "QC_histograms_allSpots.pdf"), width = 8, height = 2.5)
par(mfrow = c(1, 4))
hist(colData(spe)$sum_umi, xlab = "sum", main = "UMIs per spot all donors")
hist(colData(spe)$sum_gene, xlab = "detected", main = "Genes per spot all donors")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs all donors")
hist(colData(spe)$segmentation_info, xlab = "no. cells", main = "No. cells per spot all donors")
par(mfrow = c(1, 1))
dev.off()

# Create additional QC Metrics
qcstats <- perCellQCMetrics(spe, subsets = list(
  Mito = which(seqnames(spe) == "chrM"))
    )

qcfilter <- DataFrame(
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher")
    )

colSums(as.matrix(qcfilter))

# check thresholds
thresh_low_lib_size <- attr(qcfilter$low_lib_size, "thresholds")["lower"]
thresh_low_lib_size  #105.2499
thresh_low_n_features <- attr(qcfilter$low_n_features, "thresholds")["lower"]
thresh_low_n_features  #78.47502
thresh_high_subsets_mito_percent <- attr(qcfilter$high_subsets_Mito_percent, "thresholds")["higher"]
thresh_high_subsets_mito_percent  #45.41988

qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent

# Add QC metrics to spe object
spe$scran_low_lib_size <-
  factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
spe$scran_low_n_features <-
  factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
spe$scran_high_subsets_Mito_percent <-
  factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))
spe$scran_discard <-
  factor(qcfilter$discard, levels = c("TRUE", "FALSE"))

# QC plot of thresholds
pdf(file = here::here("plots", "QC_plots", "QC_thresholds_allSpots.pdf"), width = 8, height = 2.5)
plotQC(spe,
    type = "scatter",
    metric_x = "segmentation_info",
    metric_y = "sum_umi",
    threshold_y = thresh_low_lib_size) +
    ggtitle("Sum UMI vs. Segmentation Count all donors")
plotQC(spe,
    type = "scatter",
    metric_x = "segmentation_info",
    metric_y = "sum_gene",
    threshold_y = thresh_low_n_features) +
    ggtitle("Sum of Genes vs. Segmentation Count all donors")
plotQC(spe,
    type = "scatter",
    metric_x = "segmentation_info",
    metric_y = "subsets_mito_percent",
    threshold_y = thresh_high_subsets_mito_percent) +
    ggtitle("Mito Percent vs. Segmentation Count all donors")
dev.off()

# QC plot of tissue spots discarded
pdf(file = here::here("plots", "QC_plots", "QC_Discarded_Spots.pdf"), width = 7, height = 6.75)
df <- cbind.data.frame(colData(spe), spatialCoords(spe))
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) +
    facet_wrap(~ sample_id, nrow = 3, scales = "free") +
    geom_point(aes(color = in_tissue), size = 0.1) +
    scale_color_manual(values = "gray85") +
    new_scale_color() +
    geom_point(data = df[df$scran_discard, , drop = FALSE],
        aes(color = scran_discard), size = 0.1) +
    scale_color_manual(values = "red") +
    scale_y_reverse() +
    ggtitle("Spot-level QC") +
    theme_bw() +
    theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

save(spe, file = here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
