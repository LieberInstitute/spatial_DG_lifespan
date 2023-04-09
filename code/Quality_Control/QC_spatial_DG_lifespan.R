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
spe <-
    readRDS(here::here("processed-data", "02_build_spe", "spe.rds"))

# Check spatial orientation of Br8533 with spatialLIBD
vis_gene(
    spe = spe,
    sampleid = "Br8533",
    geneid = rownames(spe)[which(rowData(spe)$gene_name == "TNNT2")],
    assayname = "counts",
    minCount = 0,
    viridis = FALSE,
    alpha = 0.5,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br8533",
    clustervar = "10x_kmeans_3_clusters",
    minCount = 0,
    alpha = 0.5,
    point_size = 2
)

## Image transformations

spe <- rotateImg(spe, sample_id = "Br8533", degrees = 90)
spe <- mirrorImg(spe, sample_id = "Br8533", axis = "v")

# Check spatial orientation of Br8533 with spatialLIBD
vis_gene(
    spe = spe,
    sampleid = "Br8533",
    geneid = rownames(spe)[which(rowData(spe)$gene_name == "TNNT2")],
    assayname = "counts",
    minCount = 0,
    viridis = FALSE,
    alpha = 0.5,
    point_size = 2
)

vis_clus(
    spe = spe,
    sampleid = "Br8533",
    clustervar = "10x_kmeans_3_clusters",
    minCount = 0,
    alpha = 0.5,
    point_size = 2
)


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
    point_size = 2,
    alpha = 0.5
)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))

# Create QC Metrics
qcstats <- perCellQCMetrics(spe, subsets = list(Mito = which(seqnames(spe) == "chrM")))

qcfilter <- DataFrame(
    low_lib_size = isOutlier(
        qcstats$sum,
        type = "lower",
        log = TRUE,
        batch = spe$sample_id
    ),
    low_n_features = isOutlier(
        qcstats$detected,
        type = "lower",
        log = TRUE,
        batch = spe$sample_id
    ),
    high_subsets_Mito_percent = isOutlier(
        qcstats$subsets_Mito_percent,
        type = "higher",
        batch = spe$sample_id
    )
)

colSums(as.matrix(qcfilter))

qcfilter$discard <-
    qcfilter$low_lib_size | qcfilter$low_n_features |
    qcfilter$high_subsets_Mito_percent

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
pdf(
    file = here::here("plots", "QC_plots", "QC_histograms_allSpots.pdf"),
    width = 8.5,
    height = 2.5
)
par(mfrow = c(1, 4))
hist(colData(spe)$sum_umi, xlab = "sum", main = "UMIs per spot all donors")
hist(colData(spe)$sum_gene, xlab = "detected", main = "Genes per spot all donors")
hist(colData(spe)$subsets_mito_percent,
    xlab = "percent mito",
    main = "Percent mito all donors")
hist(colData(spe)$NBW, xlab = "no. cells", main = "No. cells per spot all donors")
par(mfrow = c(1, 1))
dev.off()

# Violin Plots of QC metrics
pdf(here("plots",
    "QC_plots",
    "QC_Violinplots_capture_area.pdf"),
    width = 15)

## Mito rate
plotColData(spe,
    x = "sample_id",
    y = "subsets_mito_percent",
    colour_by = "scran_high_subsets_Mito_percent") +
    ggtitle("Mito Precent by capture area")

## low sum
plotColData(spe,
    x = "sample_id",
    y = "sum",
    colour_by = "scran_low_lib_size") +
    scale_y_log10() +
    ggtitle("Total UMI count by capture area")

## low detected
plotColData(spe,
    x = "sample_id",
    y = "detected",
    colour_by = "scran_low_n_features") +
    scale_y_log10() +
    ggtitle("# detected genes by capture area")

# Mito rate vs n detected features
plotColData(
    spe,
    x = "detected",
    y = "subsets_mito_percent",
    colour_by = "scran_discard",
    point_size = 2.5,
    point_alpha = 0.5
) +
    ggtitle("Mito Precent by detected genes")

dev.off()

# QC plot of tissue spots discarded
for (i in colnames(qcfilter)) {
    vis_grid_clus(
        spe = spe,
        clustervar = paste0("scran_", i),
        pdf = here::here("plots", "QC_plots", paste0("scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        point_size = 2,
        alpha = 0.5
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

# Upload manual annotations
ManualA <- read.csv(file = here("processed-data","spatialLIBD_manual_annotations",
    "HPC annotations Lex033023.csv"))

stopifnot(ManualA$spot_name == colnames(spe))

ManualA <- as.vector(ManualA$ManualAnnotation)

spe$ManualAnnotation <- ManualA

spe$ManualAnnotation <- as.factor(spe$ManualAnnotation)

# Drop Br3874 as it has no identifiable dentate gyrus
spe <- spe[, !spe$sample_id %in% c("Br3874")]

# Create directory for QC processed spe
dir_rdata <- here::here("processed-data", "QC_processed_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

saveRDS(spe,
    file = here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()


#─ Session info ──────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.1 Patched (2022-08-30 r82775)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-03-27
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/bin/pandoc
#
#─ Packages ──────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# AnnotationDbi            1.58.0    2022-04-26 [2] Bioconductor
# AnnotationHub            3.4.0     2022-04-26 [2] Bioconductor
# assertthat               0.2.1     2019-03-21 [2] CRAN (R 4.2.1)
# attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
# beachmat                 2.12.0    2022-04-26 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
# benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
# Biobase                * 2.56.0    2022-04-26 [2] Bioconductor
# BiocFileCache            2.4.0     2022-04-26 [2] Bioconductor
# BiocGenerics           * 0.42.0    2022-04-26 [2] Bioconductor
# BiocIO                   1.6.0     2022-04-26 [2] Bioconductor
# BiocManager              1.30.19   2022-10-25 [1] CRAN (R 4.2.1)
# BiocNeighbors            1.14.0    2022-04-26 [2] Bioconductor
# BiocParallel             1.30.3    2022-06-05 [2] Bioconductor
# BiocSingular             1.12.0    2022-04-26 [2] Bioconductor
# BiocVersion              3.15.2    2022-03-29 [2] Bioconductor
# Biostrings               2.64.1    2022-08-18 [2] Bioconductor
# bit                      4.0.4     2020-08-04 [2] CRAN (R 4.2.1)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                     1.2.3     2022-04-10 [2] CRAN (R 4.2.1)
# bluster                  1.6.0     2022-04-26 [2] Bioconductor
# bslib                    0.4.0     2022-07-16 [2] CRAN (R 4.2.1)
# cachem                   1.0.6     2021-08-19 [2] CRAN (R 4.2.1)
# cli                      3.4.1     2022-09-23 [1] CRAN (R 4.2.1)
# cluster                  2.1.4     2022-08-22 [3] CRAN (R 4.2.1)
# codetools                0.2-18    2020-11-04 [3] CRAN (R 4.2.1)
# colorspace               2.0-3     2022-02-21 [2] CRAN (R 4.2.1)
# config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                   1.5.1     2022-03-26 [2] CRAN (R 4.2.1)
# curl                     4.3.3     2022-10-06 [1] CRAN (R 4.2.1)
# data.table               1.14.2    2021-09-27 [2] CRAN (R 4.2.1)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                   2.2.1     2022-06-27 [2] CRAN (R 4.2.1)
# DelayedArray             0.22.0    2022-04-26 [2] Bioconductor
# DelayedMatrixStats       1.18.0    2022-04-26 [2] Bioconductor
# desc                     1.4.1     2022-03-06 [2] CRAN (R 4.2.1)
# digest                   0.6.30    2022-10-18 [1] CRAN (R 4.2.1)
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dotCall64                1.0-2     2022-10-03 [1] CRAN (R 4.2.1)
# dplyr                    1.0.10    2022-09-01 [1] CRAN (R 4.2.1)
# dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils             1.16.0    2022-04-26 [2] Bioconductor
# DT                       0.24      2022-08-09 [2] CRAN (R 4.2.1)
# edgeR                    3.38.4    2022-08-07 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
# ExperimentHub            2.4.0     2022-04-26 [2] Bioconductor
# fansi                    1.0.3     2022-03-24 [2] CRAN (R 4.2.1)
# farver                   2.1.1     2022-07-06 [2] CRAN (R 4.2.1)
# fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.2.1)
# fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# fs                       1.5.2     2021-12-08 [2] CRAN (R 4.2.1)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           * 1.32.3    2022-08-09 [2] Bioconductor
# GenomeInfoDbData         1.2.8     2022-08-30 [2] Bioconductor
# GenomicAlignments        1.32.1    2022-07-24 [2] Bioconductor
# GenomicRanges          * 1.48.0    2022-04-26 [2] Bioconductor
# ggbeeswarm               0.6.0     2017-08-07 [2] CRAN (R 4.2.1)
# ggnewscale             * 0.4.8     2022-10-06 [1] CRAN (R 4.2.1)
# ggplot2                * 3.4.0     2022-11-04 [1] CRAN (R 4.2.1)
# ggrepel                  0.9.1     2021-01-15 [2] CRAN (R 4.2.1)
# ggside                   0.2.1     2022-07-20 [2] CRAN (R 4.2.1)
# ggspavis               * 1.2.0     2022-04-26 [2] Bioconductor
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# golem                    0.3.3     2022-07-13 [2] CRAN (R 4.2.1)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                   0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
# HDF5Array                1.24.2    2022-08-02 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# htmltools                0.5.3     2022-07-18 [2] CRAN (R 4.2.1)
# htmlwidgets              1.5.4     2021-09-08 [2] CRAN (R 4.2.1)
# httpuv                   1.6.5     2022-01-05 [2] CRAN (R 4.2.1)
# httr                     1.4.4     2022-08-17 [2] CRAN (R 4.2.1)
# igraph                   1.3.4     2022-07-19 [2] CRAN (R 4.2.1)
# interactiveDisplayBase   1.34.0    2022-04-26 [2] Bioconductor
# IRanges                * 2.30.1    2022-08-18 [2] Bioconductor
# irlba                    2.3.5     2021-12-06 [2] CRAN (R 4.2.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
# jsonlite                 1.8.3     2022-10-21 [1] CRAN (R 4.2.1)
# KEGGREST                 1.36.3    2022-07-12 [2] Bioconductor
# knitr                    1.40      2022-08-24 [2] CRAN (R 4.2.1)
# labeling                 0.4.2     2020-10-20 [2] CRAN (R 4.2.1)
# later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
# lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.1)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.2.1)
# limma                    3.52.2    2022-06-19 [2] Bioconductor
# locfit                   1.5-9.6   2022-07-11 [2] CRAN (R 4.2.1)
# magick                   2.7.3     2021-08-18 [2] CRAN (R 4.2.1)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# maps                     3.4.1     2022-10-30 [1] CRAN (R 4.2.1)
# Matrix                   1.5-1     2022-09-13 [1] CRAN (R 4.2.1)
# MatrixGenerics         * 1.8.1     2022-06-26 [2] Bioconductor
# matrixStats            * 0.62.0    2022-04-19 [2] CRAN (R 4.2.1)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# metapod                  1.4.0     2022-04-26 [2] Bioconductor
# mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# paletteer                1.4.1     2022-08-15 [2] CRAN (R 4.2.1)
# pillar                   1.8.1     2022-08-19 [2] CRAN (R 4.2.1)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# pkgload                  1.3.0     2022-06-27 [2] CRAN (R 4.2.1)
# plotly                   4.10.0    2021-10-09 [2] CRAN (R 4.2.1)
# png                      0.1-7     2013-12-03 [2] CRAN (R 4.2.1)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
# purrr                    0.3.5     2022-10-06 [1] CRAN (R 4.2.1)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                  2.12.0    2022-06-28 [2] CRAN (R 4.2.1)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                     1.0.9     2022-07-08 [2] CRAN (R 4.2.1)
# RCurl                    1.98-1.8  2022-07-30 [2] CRAN (R 4.2.1)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rhdf5                    2.40.0    2022-04-26 [2] Bioconductor
# rhdf5filters             1.8.0     2022-04-26 [2] Bioconductor
# Rhdf5lib                 1.18.2    2022-05-15 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                    1.0.6     2022-09-24 [1] CRAN (R 4.2.1)
# roxygen2                 7.2.1     2022-07-18 [2] CRAN (R 4.2.1)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools                2.12.0    2022-04-26 [2] Bioconductor
# RSQLite                  2.2.16    2022-08-17 [2] CRAN (R 4.2.1)
# rstudioapi               0.14      2022-08-22 [2] CRAN (R 4.2.1)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer              1.56.1    2022-06-23 [2] Bioconductor
# S4Vectors              * 0.34.0    2022-04-26 [2] Bioconductor
# sass                     0.4.2     2022-07-16 [2] CRAN (R 4.2.1)
# ScaledMatrix             1.4.0     2022-04-26 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater                 * 1.24.0    2022-04-26 [2] Bioconductor
# scran                  * 1.24.0    2022-04-26 [2] Bioconductor
# scuttle                * 1.6.3     2022-08-23 [2] Bioconductor
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shiny                    1.7.2     2022-07-19 [2] CRAN (R 4.2.1)
# shinyWidgets             0.7.2     2022-08-07 [2] CRAN (R 4.2.1)
# SingleCellExperiment   * 1.18.0    2022-04-26 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
# sparseMatrixStats        1.8.0     2022-04-26 [2] Bioconductor
# SpatialExperiment      * 1.6.1     2022-08-09 [2] Bioconductor
# spatialLIBD            * 1.8.10    2022-07-26 [2] Bioconductor
# statmod                  1.4.37    2022-08-12 [2] CRAN (R 4.2.1)
# stringi                  1.7.8     2022-07-11 [2] CRAN (R 4.2.1)
# stringr                  1.4.1     2022-08-20 [2] CRAN (R 4.2.1)
# SummarizedExperiment   * 1.26.1    2022-04-29 [2] Bioconductor
# tibble                   3.1.8     2022-07-22 [2] CRAN (R 4.2.1)
# tidyr                    1.2.1     2022-09-08 [1] CRAN (R 4.2.1)
# tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.2.1)
# usethis                  2.1.6     2022-05-25 [2] CRAN (R 4.2.1)
# utf8                     1.2.2     2021-07-24 [2] CRAN (R 4.2.1)
# vctrs                    0.5.0     2022-10-22 [1] CRAN (R 4.2.1)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite              0.4.1     2022-08-22 [2] CRAN (R 4.2.1)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xfun                     0.32      2022-08-10 [2] CRAN (R 4.2.1)
# XML                      3.99-0.10 2022-06-09 [2] CRAN (R 4.2.1)
# xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
# XVector                  0.36.0    2022-04-26 [2] Bioconductor
# yaml                     2.3.5     2022-02-21 [2] CRAN (R 4.2.1)
# zlibbioc                 1.42.0    2022-04-26 [2] Bioconductor
#
# [1] /users/aramnaut/R/4.2
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/R/4.2/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/R/4.2/lib64/R/library
#
#────────────────────────────────────────────────────────────────────────────────────────────────────
#
