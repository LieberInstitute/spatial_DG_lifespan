######################################
# spatial_DG_lifespan project
# SPE to Seurat conversion for PRECAST
# Anthony Ramnauth, March 30 2023
######################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
    library(Seurat)
    library(SeuratData)
    library(sessioninfo)
})

# Load SPE
spe <-
    readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

# Add row and col columns in colData(spe), specific names needed for meta.data
spe$row <- spe$array_row
spe$col <- spe$array_col

DG_seu <- CreateSeuratObject(
      counts=as.matrix(counts(spe)),
      meta.data=data.frame(colData(spe)),
      project="DG_lifespan")

seuList <- list()

samples <- unique(spe$sample_id)

for (i in seq_along(samples)){
seuList[[i]] = subset(x = DG_seu, subset = sample_id == samples[i])
}

save(seuList, file = here("processed-data", "Seurat", "seuList.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

#[1] "Reproducibility information:"
#> Sys.time()
#[1] "2023-03-30 14:45:37 EDT"
#> proc.time()
#    user   system  elapsed
# 546.059   43.380 1206.262
#> options(width = 120)
#> session_info()
#- Session info --------------------------------------------------------------------------------------
# setting  value
# version  R version 4.2.1 Patched (2022-08-30 r82775)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-03-30
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/bin/pandoc
#
#- Packages ------------------------------------------------------------------------------------------
# package              * version  date (UTC) lib source
# abind                  1.4-5    2016-07-21 [2] CRAN (R 4.2.1)
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.2.1)
# beachmat               2.12.0   2022-04-26 [2] Bioconductor
# Biobase              * 2.56.0   2022-04-26 [2] Bioconductor
# BiocGenerics         * 0.42.0   2022-04-26 [2] Bioconductor
# BiocParallel           1.30.3   2022-06-05 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.4.1    2022-09-23 [1] CRAN (R 4.2.1)
# cluster                2.1.4    2022-08-22 [3] CRAN (R 4.2.1)
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.2.1)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.2.1)
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.2.1)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.2.1)
# data.table             1.14.2   2021-09-27 [2] CRAN (R 4.2.1)
# DBI                    1.1.3    2022-06-18 [2] CRAN (R 4.2.1)
# DelayedArray           0.22.0   2022-04-26 [2] Bioconductor
# DelayedMatrixStats     1.18.0   2022-04-26 [2] Bioconductor
# deldir                 1.0-6    2021-10-23 [2] CRAN (R 4.2.1)
# digest                 0.6.30   2022-10-18 [1] CRAN (R 4.2.1)
# dplyr                  1.0.10   2022-09-01 [1] CRAN (R 4.2.1)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils           1.16.0   2022-04-26 [2] Bioconductor
# edgeR                  3.38.4   2022-08-07 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.2.1)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.2.1)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.2.1)
# fitdistrplus           1.1-8    2022-03-10 [1] CRAN (R 4.2.1)
# future                 1.29.0   2022-11-06 [1] CRAN (R 4.2.1)
# future.apply           1.10.0   2022-11-05 [1] CRAN (R 4.2.1)
# generics               0.1.3    2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.32.3   2022-08-09 [2] Bioconductor
# GenomeInfoDbData       1.2.8    2022-08-30 [2] Bioconductor
# GenomicRanges        * 1.48.0   2022-04-26 [2] Bioconductor
# ggplot2                3.4.0    2022-11-04 [1] CRAN (R 4.2.1)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.2.1)
# ggridges               0.5.3    2021-01-08 [2] CRAN (R 4.2.1)
# globals                0.16.1   2022-08-28 [2] CRAN (R 4.2.1)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.2.1)
# goftest                1.2-3    2021-10-07 [1] CRAN (R 4.2.1)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.2.1)
# gtable                 0.3.1    2022-09-01 [1] CRAN (R 4.2.1)
# HDF5Array              1.24.2   2022-08-02 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.2.1)
# htmltools              0.5.3    2022-07-18 [2] CRAN (R 4.2.1)
# htmlwidgets            1.5.4    2021-09-08 [2] CRAN (R 4.2.1)
# httpuv                 1.6.5    2022-01-05 [2] CRAN (R 4.2.1)
# httr                   1.4.4    2022-08-17 [2] CRAN (R 4.2.1)
# ica                    1.0-3    2022-07-08 [1] CRAN (R 4.2.1)
# igraph                 1.3.4    2022-07-19 [2] CRAN (R 4.2.1)
# IRanges              * 2.30.1   2022-08-18 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.2.1)
# jsonlite               1.8.3    2022-10-21 [1] CRAN (R 4.2.1)
# KernSmooth             2.23-20  2021-05-03 [3] CRAN (R 4.2.1)
# later                  1.3.0    2021-08-18 [2] CRAN (R 4.2.1)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.2.1)
# lazyeval               0.2.2    2019-03-15 [2] CRAN (R 4.2.1)
# leiden                 0.4.3    2022-09-10 [1] CRAN (R 4.2.1)
# lifecycle              1.0.3    2022-10-07 [1] CRAN (R 4.2.1)
# limma                  3.52.2   2022-06-19 [2] Bioconductor
# listenv                0.8.0    2019-12-05 [2] CRAN (R 4.2.1)
# lmtest                 0.9-40   2022-03-21 [2] CRAN (R 4.2.1)
# locfit                 1.5-9.6  2022-07-11 [2] CRAN (R 4.2.1)
# magick                 2.7.3    2021-08-18 [2] CRAN (R 4.2.1)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.2.1)
# MASS                   7.3-58.1 2022-08-03 [3] CRAN (R 4.2.1)
# Matrix                 1.5-1    2022-09-13 [1] CRAN (R 4.2.1)
# MatrixGenerics       * 1.8.1    2022-06-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.2.1)
# mime                   0.12     2021-09-28 [2] CRAN (R 4.2.1)
# miniUI                 0.1.1.1  2018-05-18 [2] CRAN (R 4.2.1)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.2.1)
# nlme                   3.1-158  2022-06-15 [3] CRAN (R 4.2.1)
# parallelly             1.32.1   2022-07-21 [2] CRAN (R 4.2.1)
# patchwork              1.1.2    2022-08-19 [2] CRAN (R 4.2.1)
# pbapply                1.5-0    2021-09-16 [2] CRAN (R 4.2.1)
# pillar                 1.8.1    2022-08-19 [2] CRAN (R 4.2.1)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.2.1)
# plotly                 4.10.0   2021-10-09 [2] CRAN (R 4.2.1)
# plyr                   1.8.7    2022-03-24 [2] CRAN (R 4.2.1)
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.2.1)
# polyclip               1.10-4   2022-10-20 [1] CRAN (R 4.2.1)
# progressr              0.11.0   2022-09-02 [1] CRAN (R 4.2.1)
# promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.2.1)
# purrr                  0.3.5    2022-10-06 [1] CRAN (R 4.2.1)
# R.methodsS3            1.8.2    2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                   1.25.0   2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                2.12.0   2022-06-28 [2] CRAN (R 4.2.1)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.2.1)
# RANN                   2.6.1    2019-01-08 [2] CRAN (R 4.2.1)
# rappdirs               0.3.3    2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.9    2022-07-08 [2] CRAN (R 4.2.1)
# RcppAnnoy              0.0.19   2021-07-30 [2] CRAN (R 4.2.1)
# RCurl                  1.98-1.8 2022-07-30 [2] CRAN (R 4.2.1)
# reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.2.1)
# reticulate             1.25     2022-05-11 [2] CRAN (R 4.2.1)
# rhdf5                  2.40.0   2022-04-26 [2] Bioconductor
# rhdf5filters           1.8.0    2022-04-26 [2] Bioconductor
# Rhdf5lib               1.18.2   2022-05-15 [2] Bioconductor
# rjson                  0.2.21   2022-01-09 [2] CRAN (R 4.2.1)
# rlang                  1.0.6    2022-09-24 [1] CRAN (R 4.2.1)
# ROCR                   1.0-11   2020-05-02 [2] CRAN (R 4.2.1)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.2.1)
# Rtsne                  0.16     2022-04-17 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.34.0   2022-04-26 [2] Bioconductor
# scales                 1.2.1    2022-08-20 [2] CRAN (R 4.2.1)
# scattermore            0.8      2022-02-14 [1] CRAN (R 4.2.1)
# sctransform            0.3.5    2022-09-21 [1] CRAN (R 4.2.1)
# scuttle                1.6.3    2022-08-23 [2] Bioconductor
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.2.1)
# Seurat               * 4.3.0    2022-11-18 [1] CRAN (R 4.2.1)
# SeuratData           * 0.2.2    2023-03-30 [1] Github (satijalab/seurat-data@d6a8ce6)
# SeuratObject         * 4.1.3    2022-11-07 [1] CRAN (R 4.2.1)
# shiny                  1.7.2    2022-07-19 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.18.0   2022-04-26 [2] Bioconductor
# sp                     1.5-1    2022-11-07 [1] CRAN (R 4.2.1)
# sparseMatrixStats      1.8.0    2022-04-26 [2] Bioconductor
# SpatialExperiment    * 1.6.1    2022-08-09 [2] Bioconductor
# spatstat.data          3.0-0    2022-10-21 [1] CRAN (R 4.2.1)
# spatstat.explore       3.0-5    2022-11-10 [1] CRAN (R 4.2.1)
# spatstat.geom          3.0-3    2022-10-25 [1] CRAN (R 4.2.1)
# spatstat.random        3.0-1    2022-11-03 [1] CRAN (R 4.2.1)
# spatstat.sparse        3.0-0    2022-10-21 [1] CRAN (R 4.2.1)
# spatstat.utils         3.0-1    2022-10-19 [1] CRAN (R 4.2.1)
# stringi                1.7.8    2022-07-11 [2] CRAN (R 4.2.1)
# stringr                1.4.1    2022-08-20 [2] CRAN (R 4.2.1)
# SummarizedExperiment * 1.26.1   2022-04-29 [2] Bioconductor
# survival               3.4-0    2022-08-09 [3] CRAN (R 4.2.1)
# tensor                 1.5      2012-05-05 [1] CRAN (R 4.2.1)
# tibble                 3.1.8    2022-07-22 [2] CRAN (R 4.2.1)
# tidyr                  1.2.1    2022-09-08 [1] CRAN (R 4.2.1)
# tidyselect             1.2.0    2022-10-10 [1] CRAN (R 4.2.1)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.2.1)
# uwot                   0.1.14   2022-08-22 [2] CRAN (R 4.2.1)
# vctrs                  0.5.0    2022-10-22 [1] CRAN (R 4.2.1)
# viridisLite            0.4.1    2022-08-22 [2] CRAN (R 4.2.1)
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.2.1)
# XVector                0.36.0   2022-04-26 [2] Bioconductor
# zlibbioc               1.42.0   2022-04-26 [2] Bioconductor
# zoo                    1.8-10   2022-04-15 [2] CRAN (R 4.2.1)
#
# [1] /users/aramnaut/R/4.2
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/R/4.2/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/R/4.2/lib64/R/library
#
#-----------------------------------------------------------------------------------------------------
