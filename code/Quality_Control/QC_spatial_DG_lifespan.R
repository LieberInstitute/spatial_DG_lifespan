##############################################################
# spatial_DG_lifespan project
# Script for quality control (QC) & Normalization (log counts)
# Anthony Ramnauth, Apr 18 2022
##############################################################

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
    library(sessioninfo)
})

# load saved SPE object
spe <- readRDS(here::here("processed-data", "02_build_spe", "spe.rds"))

# Create directory for QC plots
dir_rdata <- here::here("plots", "QC_plots")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

dim(spe)

# check SPE object contains only spots over tissue
table(colData(spe)$in_tissue)
all(colData(spe)$in_tissue)


# Plot spatial coordinates & orientation of within tissue spots
vis_grid_clus(
    spe = spe,
    clustervar = "in_tissue",
    pdf = here::here("plots", "QC_plots", "in_tissue_grid.pdf"),
    sort_clust = FALSE,
    colors = c("TRUE" = "grey90", "FALSE" = "orange"),
    point_size = 2
)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))

# Create QC Metrics
qcstats <- perCellQCMetrics(spe, subsets = list(
    Mito = which(seqnames(spe) == "chrM")
))

qcfilter <- DataFrame(
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = spe$sample_id),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = spe$sample_id),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher", batch = spe$sample_id)
)

colSums(as.matrix(qcfilter))

qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent

# Add QC metrics to spe object
spe$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
spe$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
spe$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
spe$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

# plot histograms of QC metrics
pdf(file = here::here("plots", "QC_plots", "QC_histograms_allSpots.pdf"), width = 8.5, height = 2.5)
par(mfrow = c(1, 4))
hist(colData(spe)$sum_umi, xlab = "sum", main = "UMIs per spot all donors")
hist(colData(spe)$sum_gene, xlab = "detected", main = "Genes per spot all donors")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito all donors")
hist(colData(spe)$count, xlab = "no. cells", main = "No. cells per spot all donors")
par(mfrow = c(1, 1))
dev.off()

# QC plot of tissue spots discarded
for (i in colnames(qcfilter)) {
    vis_grid_clus(
        spe = spe,
        clustervar = paste0("scran_", i),
        pdf = here::here("plots", "QC_plots", paste0("scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        point_size = 2
    )
}

# remove combined set of low-quality spots

spe <- spe[, colData(spe)$scran_discard == FALSE]
dim(spe)

# calculate logcounts (log-transformed normalized counts) and store in object
spe$scran_quick_cluster <- quickCluster(spe, block = spe$sample_id)
spe <- computeSumFactors(spe, clusters = spe$scran_quick_cluster)
table(spe$scran_quick_cluster)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)

spe <- logNormCounts(spe)

# check
assayNames(spe)

# Create directory for QC processed spe
dir_rdata <- here::here("processed-data", "QC_processed_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

saveRDS(spe, file = here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
