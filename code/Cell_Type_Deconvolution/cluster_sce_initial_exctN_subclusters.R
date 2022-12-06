#####################################################
# spatial_DG_lifespan project
# Sub-clustering of sce object for excitatory neurons
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

# Subset for Excitatory Neurons

sce_full <- sce

# select excitatory neuron clusters
# (identified based on marker expression heatmap from previous script)
clus_select <- c("ExctN")

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

# secondary clustering of excitatory neurons

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and graph-based clustering

set.seed(12345)
clus <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 2000),
    second = NNGraphParam(k = 25)
  )
)

table(clus)
#clus
#    1     2     3     4     5     6     7
#13228  8517 24607  1429  5060  4465  5637

colLabels(sce) <- clus

table(colLabels(sce), colData(sce)$Dataset)
#    Franjic_etal_2022 Zhong_etal_2020 Zhou_etal_2022
#  1              3175               1          10052
#  2                88            8418             11
#  3             10491              75          14041
#  4                23               0           1406
#  5                 4               0           5056
#  6               201               0           4264
#  7              5453             183              1

# Check marker genes violin plots

pdf(file = here::here("plots", "sce_plots", "ExctN_cluster_markers_sce.pdf"))

plotExpression(sce, features=c("PROX1", "CALB1", "PDLIM5", "SGCZ", "ARHGAP24",
    "DLC1", "CFAP299", "SYN3", "HGF", "ACVR1C", "SYT13", "ROBO1", "COL5A2"),
    x="label", colour_by="label")

dev.off()

# Plot UMAP after clustering
pdf(file = here::here("plots", "sce_plots", "DG_UMAP_ExctN_initial_clusters_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label",
    point_alpha = 0.3, point_size = 0.5) +
  ggtitle("Unsupervised clustering of initial ExctN clusters")

dev.off()

# Create heatmap for markers of mean expresion with z-scores
markers <- c(
    # GC
    "PROX1", "CALB1", "PDLIM5", "SGCZ",
    # im GC
    "DCX", "BHLHE22", "STMN1",
    # Mossy Cells
    "ARHGAP24", "DLC1", "ADCYAP1",
    # CA3 PNs
    "CFAP299", "SYN3",
    # CA2 PNs
    "HGF",
    # CA1 PNs
    "ACVR1C", "SYT13", "GRIK1", "ACVR1C",
    # Sub PNs
    "ROBO1", "COL5A2", "FN1",
    # Progenitors
    "PAX6", "HOPX", "EOMES", "SOX2", "NEUROG1")

# marker labels
marker_labels <- c(
  rep("Granular_cells", 4),
  rep("IM_Granular_cells", 3),
  rep("Mossy_Cells", 3),
  rep("CA3", 2),
  rep("CA2", 1),
  rep("CA1", 4),
  rep("Sub", 3),
  rep("Progen", 5))

marker_labels <-
  factor(marker_labels, levels = unique(marker_labels))

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

colors_markers <- list(marker = c(
  Granular_cells = "deepskyblue",
  IM_Granular_cells = "dodgerblue",
  Mossy_Cells = "blue",
  CA3 = "deepskyblue3",
  CA2 = "blue3",
  CA1 = "deepskyblue4",
  Sub = "darkblue",
  Progen = "cyan"))

# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels,
  show_annotation_name = FALSE,
  show_legend = TRUE,
  col = colors_markers
  )

pdf(file = here::here("plots", "sce_plots", "Initial_Heatmap_ExctN_markers_sce.pdf"), width = 12, height = 8)

Heatmap(
  hm_mat,
  name = "z-score",
  column_title = "ExctN markers",
  column_title_gp = gpar(fontface = "bold"),
  bottom_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_split = marker_labels,
  column_names_gp = gpar(fontface = "italic"),
  rect_gp = gpar(col = "gray50", lwd = 0.5)
    )

dev.off()

# 1 = Sub, 2 = Progen, 3 = Mixed, 4 = Mixed, 5 = Mixed, 6 = Mixed, 7 = Mixed
# Sub & Progen are the only obvious clusters

######################
# Store cluster labels
######################

# store secondary clustering labels in full SCE object

# check unique barcode IDs
table(duplicated(colnames(sce)))
table(duplicated(colnames(sce_full)))
table(colnames(sce) %in% colnames(sce_full))

# match and store cluster labels
clus_excitatory <- rep(NA, ncol(sce_full))
names(clus_excitatory) <- colnames(sce_full)
clus_excitatory[colnames(sce)] <- colData(sce)$label

colData(sce_full)$label_excitatory <- clus_excitatory

# check
table(colData(sce_full)$label)
table(colData(sce_full)$label_excitatory)
table(colData(sce_full)$label_excitatory, useNA = "always")

# Add character to label_Mixed to denote cluster number & from Mixed

colData(sce_full)$label_excitatory <-
    sub("^", "ExctN", colData(sce_full)$label_excitatory)

sce_full$label_merged[sce_full$label_merged == "ExctN"] <- NA

dfc <- data.frame(
    label_merged = colData(sce_full)$label_merged,
    label_excitatory = colData(sce_full)$label_excitatory)

dfc <- dfc %>%
    mutate(
        cell_type = coalesce(
            label_merged,
            label_excitatory))

dfc <- dfc %>%
    mutate(cell_type = recode(cell_type, ExctN1 = "Sub", ExctN2 = "NPC",
        ExctN3 =  "Mixed", ExctN4 = "Mixed", ExctN5 = "Mixed", ExctN6 = "Mixed",
        ExctN7 = "Mixed"))

sce_full$label_merged <- dfc$cell_type

# Save new sce object
saveRDS(sce_full, file = here::here("processed-data", "sce", "sce_clustered.rds"))
