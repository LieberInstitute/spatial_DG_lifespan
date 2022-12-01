###############################
# spatial_DG_lifespan project
# QC summed sce object
# Anthony Ramnauth, Nov 05 2022
###############################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scuttle)
    library(DropletUtils)
    library(here)
    library(scater)
    library(scran)
    library(scDblFinder)
    library(sessioninfo)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "build_sce", "sce_sum.rds"))

dir_plots <- here("sce_plots")

# ---------------
# Remove doublets
# ---------------

# identify and remove doublets using scDblFinder

# note: no random seed required

sce <- scDblFinder(sce, samples = "Dataset")

# number of doublets per sample
table(colData(sce)$Dataset, colData(sce)$scDblFinder.class)

# remove doublets

dim(sce)

# Violin Plots of doublets AFTER QC metrics to have y-axis be sum (UMIs)
pdf(
    file = here::here(
        "sce_plots",
        "doublets_sce_sum.pdf"
    )
)

plotColData(sce,
    x = "Dataset",
    y = "sum",
    colour_by = "scDblFinder.class") +
    scale_y_log10() +
    ggtitle("doublets of cells")

dev.off()

ix_dbl <- colData(sce)$scDblFinder.class == "doublet"
table(ix_dbl)

sce <- sce[, !ix_dbl]

dim(sce)

# ----------
# QC metrics
# ----------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$SYMBOL)
table(is_mito)
rowData(sce)$SYMBOL[is_mito]

# calculate per-spot QC metrics and store in colData
sce <- addPerCellQC(sce, subsets = list(mito = is_mito))
head(colData(sce))

colData(sce)$subsets_mito_percent[is.na(colData(sce)$subsets_mito_percent)] = 0
colData(sce)$sum[is.na(colData(sce)$sum)] = 0
colData(sce)$detected[is.na(colData(sce)$detected)] = 0

# Create QC Metrics
qcstats <- perCellQCMetrics(sce, subsets = list(mito = is_mito))

qcstats$subsets_mito_percent[is.na(qcstats$subsets_mito_percent)] = 0
qcstats$sum[is.na(qcstats$sum)] = 0
qcstats$detected[is.na(qcstats$detected)] = 0

qcfilter <- DataFrame(
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = sce$Dataset),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = sce$Dataset),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_mito_percent, type = "higher", batch = sce$Dataset)
)

colSums(as.matrix(qcfilter))

qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent

# Add QC metrics to sce object
sce$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
sce$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
sce$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
sce$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))


# plot histograms of QC metrics
pdf(
    file = here::here(
        "sce_plots",
        "histograms_sce_sum.pdf"
    ),
    width = 5,
    height = 2
)
par(mfrow = c(1, 3))
hist(colData(sce)$sum,
    xlab = "sum",
    main = "UMIs per cell")
hist(colData(sce)$detected,
    xlab = "features",
    main = "Genes per cell")
hist(colData(sce)$subsets_mito_percent,
    xlab = "percent mito",
    main = "Percent mito")
par(mfrow = c(1, 1))

dev.off()

# Violin Plots of QC metrics
pdf(
    file = here::here(
        "sce_plots",
        "violinplots_sce_sum.pdf"
    )
)
## Mito rate
plotColData(sce,
    x = "Dataset",
    y = "subsets_mito_percent",
    colour_by = "scran_high_subsets_Mito_percent") +
    ggtitle("Mito Precent by cell")

## low sum
plotColData(sce,
    x = "Dataset",
    y = "sum",
    colour_by = "scran_low_lib_size") +
    scale_y_log10() +
    ggtitle("Total UMI count by cell")

## low detected
plotColData(sce,
    x = "Dataset",
    y = "detected",
    colour_by = "scran_low_n_features") +
    scale_y_log10() +
    ggtitle("# detected genes by cell")

# Mito rate vs n detected features
plotColData(
    sce,
    x = "detected",
    y = "subsets_mito_percent",
    colour_by = "scran_discard",
    point_size = 2.5,
    point_alpha = 0.5
) +
    ggtitle("Mito Precent by detected genes")

dev.off()

# remove combined set of low-quality cells

sce <- sce[, colData(sce)$scran_discard == FALSE]

dim(sce)

# calculate logcounts (log-transformed normalized counts) and store in object
sce$scran_quick_cluster <- quickCluster(sce, block = sce$Dataset)
sce <- computeSumFactors(sce, clusters = sce$scran_quick_cluster)
table(sce$scran_quick_cluster)
summary(sizeFactors(sce))
hist(sizeFactors(sce), breaks = 40)

sce <- logNormCounts(sce)

# check
assayNames(sce)

saveRDS(sce, file = here::here("sce_objects", "QCed_sce_sum.rds"))


## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
