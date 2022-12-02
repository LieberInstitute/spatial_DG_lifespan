#####################################################
# spatial_DG_lifespan project
# Sub-clustering of sce object for inhibitory neurons
# Anthony Ramnauth, Nov 28 2022
#####################################################

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

# Subset for Inhibitory Neurons

sce_full <- sce

# select inhibitory neuron clusters
# (identified based on marker expression heatmap from previous script)
clus_select <- c("InhbN")

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

# secondary clustering of inhibitory neurons

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and Leiden clustering

set.seed(12345)
clus <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 2000),
    second = NNGraphParam(k = 20, cluster.fun = "leiden")
  )
)

table(clus)
#clus
#    1     2     3     4     5     6     7     8     9
#10500   412  2625  3269  1768  3307  1528   647    34

colLabels(sce) <- clus

table(colLabels(sce), colData(sce)$Dataset)
#    Franjic_etal_2022 Zhong_etal_2020 Zhou_etal_2022
#  1               272            5790           4438
#  2                 0             412              0
#  3               960             169           1496
#  4                 8            3258              3
#  5              1449              24            295
#  6               903             221           2183
#  7                20               0           1508
#  8                 0             647              0
#  9                34               0              0

# Check marker genes violin plots

pdf(file = here::here("plots", "sce_plots", "InhbN_cluster_markers_sce.pdf"))

plotExpression(sce, features=c("GAD1", "GAD2", "SST", "KIT", "CALB1", "CALB2", "TAC1",
    "CNR1", "PVALB", "CORT", "VIP", "NPY", "CRHBP", "CCK", "HTR3A", "NR2F2", "LAMP5", "MEIS2"),
    x="label", colour_by="label")

dev.off()

# Plot UMAP after clustering
pdf(file = here::here("plots", "sce_plots", "DG_UMAP_InhbN_clusters_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label",
    point_alpha = 0.3, point_size = 0.5) +
  ggtitle("Unsupervised clustering of InhbN clusters")

dev.off()

# Create heatmap for markers of mean expresion with z-scores
markers <- c("GAD1", "GAD2", "SST", "KIT", "CALB1", "CALB2", "TAC1","CNR1",
    "PVALB", "CORT", "VIP", "NPY", "CRHBP", "CCK", "HTR3A", "NR2F2", "LAMP5",
    "MEIS2", "DCX", "BHLHE22", "STMN1", "PAX6", "EOMES")

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

pdf(file = here::here("plots", "sce_plots", "Heatmap_InhbN_markers_sce.pdf"), width = 12, height = 8)

Heatmap(
  hm_mat,
  name = "z-score",
  column_title = "InhbN markers",
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
clus_inhibitory <- rep(NA, ncol(sce_full))
names(clus_inhibitory) <- colnames(sce_full)
clus_inhibitory[colnames(sce)] <- colData(sce)$label

colData(sce_full)$label_inhibitory <- clus_inhibitory

# check
table(colData(sce_full)$label)
table(colData(sce_full)$label_inhibitory)
table(colData(sce_full)$label_inhibitory, useNA = "always")


# Save new sce object
saveRDS(sce_full, file = here::here("processed-data", "sce", "sce_clustered.rds"))
