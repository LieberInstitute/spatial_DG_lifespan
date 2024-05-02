#######################################
# spatial_DG_lifespan project
# QC sce object from mouse neurogenesis
# Anthony Ramnauth, Jan 27 2024
#######################################

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

sce <- readRDS(here::here("sce_objects", "sce_hochgerner.rds"))

dir_plots <- here::here("sce_plots_hochgerner")

# ---------------
# Remove doublets
# ---------------

# identify and remove doublets using scDblFinder

# note: no random seed required

sce <- scDblFinder(sce, clusters = "cellType", samples = "strain")

# number of doublets per sample
table(colData(sce)$cellType, colData(sce)$scDblFinder.class)

#                     singlet doublet
#  Astro-adult             1203      29
#  Astro-juv                798      23
#  CA3-Pyr                  454      78
#  Cajal-Retzius            457      78
#  Endothelial              508      35
#  Ependymal                157      25
#  GABA                     171      38
#  GC-adult                2595      18
#  GC-juv                  3346      74
#  Immature-Astro           602      49
#  Immature-GABA            950      74
#  Immature-GC             2230     189
#  Immature-Pyr            4152     368
#  MiCajal-Retziusoglia     401      27
#  MOL                      648      56
#  Neuroblast              1321      60
#  NFOL                     220      12
#  nIPC                     319      34
#  nIPC-perin               405      63
#  OPC                      747      47
#  PVM                       75       4
#  RGL                      186      13
#  RGL_young                630      56
#  VLMC                     133      27

# remove doublets

dim(sce)

# Violin Plots of doublets AFTER QC metrics to have y-axis be sum (UMIs)
pdf(
    file = here::here("sce_plots_hochgerner",
        "doublets_sce_hochgerner.pdf"
    )
)

plotColData(sce,
    x = "cellType",
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
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = sce$strain),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = sce$strain),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_mito_percent, type = "higher", batch = sce$strain)
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
    file = here::here("sce_plots_hochgerner",
        "histograms_sce_hochgerner.pdf"
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
    file = here::here("sce_plots_hochgerner",
        "violinplots_sce_hochgerner.pdf"
    )
)
## Mito rate
plotColData(sce,
    x = "cellType",
    y = "subsets_mito_percent",
    colour_by = "scran_high_subsets_Mito_percent") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Mito Precent by cell")

## low sum
plotColData(sce,
    x = "cellType",
    y = "sum",
    colour_by = "scran_low_lib_size") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_y_log10() +
    ggtitle("Total UMI count by cell")

## low detected
plotColData(sce,
    x = "cellType",
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
sce$scran_quick_cluster <- quickCluster(sce, block = sce$strain)
sce <- computeSumFactors(sce, clusters = sce$scran_quick_cluster)
table(sce$scran_quick_cluster)
summary(sizeFactors(sce))

sce <- logNormCounts(sce)

# check
assayNames(sce)

saveRDS(sce, file = here::here("sce_objects", "QCed_sce_hochgerner.rds"))


## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
