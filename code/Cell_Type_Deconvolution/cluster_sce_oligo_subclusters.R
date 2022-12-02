###################################################
# spatial_DG_lifespan project
# Sub-clustering of sce object for oligodendrocytes
# Anthony Ramnauth, Dec 01 2022
###################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(bluster)
    library(ggplot2)
    library(viridisLite)
    library(dplyr)
    library(tidyr)
    library(harmony)
    library(ComplexHeatmap)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_clustered.rds"))

# Subset for oligodendrocytes

sce_full <- sce

# select oligo clusters
# (identified based on marker expression heatmap from previous script)
clus_select <- c("Oligo")

ix_select <- sce$label_merged %in% clus_select
table(ix_select)

sce <- sce[, ix_select]

dim(sce)

# re-calculate logcounts and PCA on subset to ensure PCs capture the most relevant variation

sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

dec <- modelGeneVar(sce, block = sce$Dataset)
top_hvgs <- getTopHVGs(dec, prop = 0.1)

set.seed(12345)
sce <- runPCA(sce, subset_row = top_hvgs, ncomponents = 50)

pca_matrix <- reducedDim(sce, "PCA")
Dataset <- colData(sce)$Dataset
stopifnot(nrow(pca_matrix) == length(Dataset))

set.seed(123)
harmony_embeddings <- HarmonyMatrix(
    pca_matrix,
    meta_data = Dataset,
    do_pca = FALSE
)

# calculate UMAP (on harmony embeddings)
set.seed(123)
sce <- runUMAP(sce, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

# secondary clustering of oligos

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and Leiden clustering

set.seed(12345)
clus <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 1000),
    second = NNGraphParam(k = 20, cluster.fun = "leiden")
  )
)

table(clus)
#clus
#    1     2     3     4
#47255 13154  3673  5686

colLabels(sce) <- clus

table(colLabels(sce), colData(sce)$Dataset)
#    Franjic_etal_2022 Zhong_etal_2020 Zhou_etal_2022
#  1             17794             970          28491
#  2             12302             779             73
#  3              3515               0            158
#  4              5281              51            354

# Check marker genes violin plots

pdf(file = here::here("plots", "sce_plots", "Oligo_cluster_markers_sce.pdf"))

plotExpression(sce, features=c("MOBP", "OLIG1", "OLIG2", "SOX10", "MOG", "CNP",
    "CLDN11", "SOX2", "PAX6", "HOPX", "EOMES"),
    x="label", colour_by="label")

dev.off()

# Plot UMAP after clustering
pdf(file = here::here("plots", "sce_plots", "DG_UMAP_Oligo_clusters_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label",
    point_alpha = 0.3, point_size = 0.5) +
  ggtitle("Unsupervised clustering of Oligo clusters")

dev.off()

# Create heatmap for markers of mean expresion with z-scores
markers <- c("MOBP", "OLIG1", "OLIG2", "SOX10", "MOG", "CNP",
    "CLDN11", "SOX2", "PAX6", "HOPX", "EOMES")

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$label)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$SYMBOL

hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers, i]))))

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

hm_mat <- t(scale_rows(t(hm_mat)))

pdf(file = here::here("plots", "sce_plots", "Heatmap_Oligo_markers_sce.pdf"), width = 12, height = 8)

Heatmap(
  hm_mat,
  name = "z-score",
  column_title = "Oligo markers",
  column_title_gp = gpar(fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_names_gp = gpar(fontface = "italic"),
  rect_gp = gpar(col = "gray50", lwd = 0.5)
    )

dev.off()

######################
# Store cluster labels
######################

# store secondary clustering labels in full SCE object

# check unique barcode IDs
table(duplicated(colnames(sce)))
table(duplicated(colnames(sce_full)))
table(colnames(sce) %in% colnames(sce_full))

# match and store cluster labels
clus_oligo <- rep(NA, ncol(sce_full))
names(clus_oligo) <- colnames(sce_full)
clus_oligo[colnames(sce)] <- colData(sce)$label

colData(sce_full)$label_oligo <- clus_oligo

# check
table(colData(sce_full)$label)
table(colData(sce_full)$label_oligo)
table(colData(sce_full)$label_oligo, useNA = "always")


# Save new sce object
saveRDS(sce_full, file = here::here("processed-data", "sce", "sce_clustered.rds"))
