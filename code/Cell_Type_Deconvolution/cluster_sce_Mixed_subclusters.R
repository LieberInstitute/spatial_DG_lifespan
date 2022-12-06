###############################################
# spatial_DG_lifespan project
# Sub-clustering of sce object for Mixed nuclei
# Anthony Ramnauth, Dec 01 2022
###############################################

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
    library(forcats)
    library(ComplexHeatmap)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_clustered.rds"))

# Subset for Mixed cells

sce_full <- sce

# select Mixed clusters
# (identified based on marker expression heatmap from previous script)
clus_select <- c("Mixed")

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

# secondary clustering of Mixed cells

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and Leiden clustering

set.seed(12345)
clus <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 1000),
    second = NNGraphParam(k = 25)
  )
)

table(clus)
#clus
#    1     2     3     4     5     6     7     8
#  305  1645   762  3510   991 10112   470  1051

colLabels(sce) <- clus

table(colLabels(sce), colData(sce)$Dataset)
#    Franjic_etal_2022 Zhong_etal_2020 Zhou_etal_2022
#  1                24               0            281
#  2              1271               2            372
#  3               460               0            302
#  4                15            3467             28
#  5                 0               0            991
#  6               561            5638           3913
#  7               269              10            191
#  8                 3               0           1048

# Check marker genes violin plots

# Plot UMAP after clustering
pdf(file = here::here("plots", "sce_plots", "DG_UMAP_initial_mixed_clusters_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label",
    point_alpha = 0.3, point_size = 0.5) +
  ggtitle("Unsupervised clustering of intiial mixed clusters")

dev.off()

# ------------
# Marker genes
# ------------

# Markers chosen from the publications of each dataset

markers <- c(
    ## Neurons
    "RBFOX3", "SNAP25", "SYT1",
    ## Excitatory Neurons
    "SLC17A7",
    # GC
    "PROX1",
    # im GC
    "DCX", "BHLHE22", "STMN1",
    # Mossy Cells
    "ARHGAP24", "DLC1",
    # CA3 PNs
    "CFAP299", "SYN3",
    # CA2 PNs
    "HGF",
    # CA1 PNs
    "ACVR1C", "SYT13",
    # Sub PNs
    "ROBO1", "COL5A2",
    # Cajal–Retzius
    "RELN",
    ## Inhibitory Neurons
    "GAD1", "GAD2",
    # inhibitory subpopulations (Some from Lukas LC & Keri Martinowich 2022-07-22)
    "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", "PVALB", "CORT", "VIP", "NPY",
    "CRHBP", "CCK", "HTR3A", "NR2F2", "LAMP5",
    # Astrocytes
    "AQP4", "GFAP", "CHRDL1",
    # Oligodendrocytes
    "MOBP",
    # macrophages / microglia
    "CD163", "C3", "PTPRC", "C1QB",
    # OPCs
    "PDGFRA", "VCAN",
    # COP
    "GPR17", "ADAM33",
    # endothelial / mural (RBPMS)
    "CLDN5", "FLT1", "RBPMS",
    # T cells
    "SKAP1", "CD247",
    # Progenitors
    "PAX6", "HOPX", "EOMES"
)

# marker labels
marker_labels <- c(
  rep("Neuron", 3),
  rep("Excitatory", 1),
  rep("Granular_cells", 1),
  rep("IM_Granular_cells", 3),
  rep("Mossy_Cells", 2),
  rep("CA3", 2),
  rep("CA2", 1),
  rep("CA1", 2),
  rep("Sub", 2),
  rep("Cajal–Retzius", 1),
  rep("Inhibitory", 17),
  rep("Astrocytes", 3),
  rep("Oligodendrocytes", 1),
  rep("Macrophages_Microglia", 4),
  rep("OPCs", 2),
  rep("COP", 2),
  rep("Endothelial_Mural", 3),
  rep("T_cells", 2),
  rep("Progenitors", 3)
)

# colors: selected from tableau20 and tableau10medium
colors_markers <- list(marker = c(
  Neuron = "blue",
  Excitatory = "darkblue",
  Granular_cells = "blue1",
  IM_Granular_cells = "blue2",
  Mossy_Cells = "blue4",
  CA3 = "dodgerblue",
  CA2 = "dodgerblue3",
  CA1 = "blue3",
  Sub = "blue4",
  'Cajal–Retzius' = "gray47",
  Inhibitory = "green",
  Astrocytes = "yellow",
  Oligodendrocytes = "plum3",
  Macrophages_Microglia = "tan",
  OPCs = "goldenrod",
  COP = "goldenrod4",
  Endothelial_Mural = "red3",
  T_cells = "tan3",
  Progenitors = "cyan"
    ))

marker_labels <-
  factor(marker_labels, levels = unique(marker_labels))

# number of nuclei per cluster
n <- table(colLabels(sce))

# heatmap data

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

# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels,
  show_annotation_name = FALSE,
  show_legend = TRUE,
  col = colors_markers
  )

pdf(file = here::here("plots", "sce_plots", "Mixed_heatmap_markers_sce.pdf"), width = 12, height = 8)

Heatmap(
  hm_mat,
  name = "z-score",
  column_title = "Mixed DG clusters",
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

# Mixed cluster
# 1 = ExctN?, 2 = Endo_M, 3 = Endo_M, 4 = InhbN, 5 = ExctN?, 6 = InhbN,
# 7 = Immune, 8 = ExctN

######################
# Store cluster labels
######################

# store secondary clustering labels in full SCE object

# check unique barcode IDs
table(duplicated(colnames(sce)))
table(duplicated(colnames(sce_full)))
table(colnames(sce) %in% colnames(sce_full))

# match and store cluster labels
clus_Mixed <- rep(NA, ncol(sce_full))
names(clus_Mixed) <- colnames(sce_full)
clus_Mixed[colnames(sce)] <- colData(sce)$label

colData(sce_full)$label_Mixed <- clus_Mixed

# check
table(colData(sce_full)$label)
table(colData(sce_full)$label_Mixed)
table(colData(sce_full)$label_Mixed, useNA = "always")

# Rename values in label_merged based on updated mixed clusters

# Add character to label_Mixed to denote cluster number & from Mixed

colData(sce_full)$label_Mixed <-
    sub("^", "Mixed", colData(sce_full)$label_Mixed)

sce_full$label_merged[sce_full$label_merged == "Mixed"] <- NA

dfc <- data.frame(
    label_merged = colData(sce_full)$label_merged,
    label_Mixed = colData(sce_full)$label_Mixed)

dfc <- dfc %>%
    mutate(
        cell_type = coalesce(
            label_merged,
            label_Mixed))

# Rename the mixed cluster values to the merged cluster values

dfc <- dfc %>%
    mutate(cell_type = recode(cell_type, Mixed1 = "ExctN", Mixed2 = "Endo_M",
        Mixed3 =  "Endo_M", Mixed4 = "InhbN", Mixed5 = "ExctN", Mixed6 = "InhbN",
        Mixed7 = "Immune", Mixed8 = "ExctN"))

sce_full$label_merged <- dfc$cell_type

# Updated UMAP

colors_clusters <- list(population = c(
  ExctN = "blue",
  InhbN = "green",
  Glia = "yellow",
  Oligo = "plum3",
  Immune = "tan",
  OPCs = "goldenrod",
  Endo_M = "red3")
    )

pdf(file = here::here("plots", "sce_plots", "Updated_merged__initial_mixed_cluster_plot_sce.pdf"))

plotReducedDim(sce_full, dimred = "UMAP.HARMONY", colour_by = "label_merged") +
  scale_color_manual(values = colors_clusters[[1]], name = "clusters (merged)") +
  theme_classic() +
  ggtitle("Combined HPC snRNAseq datasets of updated merged clustering")

dev.off()

# Save new sce object
saveRDS(sce_full, file = here::here("processed-data", "sce", "sce_clustered.rds"))
