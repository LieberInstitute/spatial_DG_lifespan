###############################
# spatial_DG_lifespan project
# QC sce object from Macaque
# Anthony Ramnauth, Dec 07 2023
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

sce <- readRDS(here::here("sce_objects", "sce_hao.rds"))

dir_plots <- here::here("sce_plots_hao")

# ---------------
# Remove doublets
# ---------------

# identify and remove doublets using scDblFinder

# note: no random seed required

sce <- scDblFinder(sce, samples = "sample_ID")

# number of doublets per sample
table(colData(sce)$sample_ID, colData(sce)$scDblFinder.class)

#            singlet doublet
#  0118_H1      6334     147
#  0118_H2      6507     289
#  0118_H3      7012     187
#  0118_H4      6327     174
#  0125_H1      3878      80
#  0614_H1      3627      74
#  0614_H2      6440     187
#  0614_H3      5675     150
#  0614_H4      4212     148
#  0625_H3      7656     355
#  0731_H1     11141     570
#  0731_H3      8643     291
#  0731_H4      9795     427
#  0801_H1-1   10641     485
#  0801_H1-2    9425     405
#  0801_H2      6640     212
#  0812_H4-1    9897     394
#  0812_H4-2    9549     360
#  1211_H2-2    9488     412
#  1211_H3-1    9244     541
#  1211_H3-2    8814     646
#  1211_H4-1    9674     584
#  1211_H4-3   10772     950
#  1211_H6-1    7924     347
#  1211_H6-2    9453     602

# remove doublets

dim(sce)

# Violin Plots of doublets AFTER QC metrics to have y-axis be sum (UMIs)
pdf(
    file = here::here("sce_plots_hao",
        "doublets_sce_hao.pdf"
    )
)

plotColData(sce,
    x = "sample_ID",
    y = "scDblFinder.score",
    colour_by = "scDblFinder.class") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90)) +
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
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)
table(is_mito)
rowData(sce)$gene_name[is_mito]

# Looks like they already removed all mitochondrial genes!

# calculate per-spot QC metrics and store in colData
sce <- addPerCellQC(sce, subsets = list(mito = is_mito))
head(colData(sce))

# Create QC Metrics
qcstats <- perCellQCMetrics(sce, subsets = list(mito = is_mito))

qcfilter <- DataFrame(
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = sce$sample_ID),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = sce$sample_ID),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_mito_percent, type = "higher", batch = sce$sample_ID)
)

colSums(as.matrix(qcfilter))

qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features | qcfilter$high_subsets_Mito_percent)

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
    file = here::here("sce_plots_hao",
        "histograms_sce_hao.pdf"
    ),
    width = 5,
    height = 2
)
par(mfrow = c(1, 3))
hist(colData(sce)$sum,
    xlab = "sum",
    main = "UMIs per cell")
hist(colData(sce)$detected,
    xlab = "detected",
    main = "Genes per cell")
hist(colData(sce)$subsets_mito_percent,
    xlab = "percent mito",
    main = "Percent mito")
par(mfrow = c(1, 1))

dev.off()

# Violin Plots of QC metrics
pdf(
    file = here::here("sce_plots_hao",
        "violinplots_sce_hao.pdf"
    )
)
## Mito rate
plotColData(sce,
    x = "sample_ID",
    y = "subsets_mito_percent",
    colour_by = "scran_high_subsets_Mito_percent") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Mito Precent by cell")

## low sum
plotColData(sce,
    x = "sample_ID",
    y = "sum",
    colour_by = "scran_low_lib_size") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_y_log10() +
    ggtitle("Total UMI count by cell")

## low detected
plotColData(sce,
    x = "sample_ID",
    y = "detected",
    colour_by = "scran_low_n_features") +
    theme(axis.text.x = element_text(angle = 90)) +
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
sce$scran_quick_cluster <- quickCluster(sce, block = sce$sample_ID)
sce <- computeSumFactors(sce, clusters = sce$scran_quick_cluster)
table(sce$scran_quick_cluster)
summary(sizeFactors(sce))

sce <- logNormCounts(sce)

# check
assayNames(sce)

saveRDS(sce, file = here::here("sce_objects", "QCed_sce_hao.rds"))


## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
