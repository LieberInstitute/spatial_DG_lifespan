######################################################
# spatial_DG_lifespan project
# Identifying clusters of sce object with marker genes
# Anthony Ramnauth, Nov 06 2022
######################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(dplyr)
    library(tidyr)
    library(forcats)
    library(ggplot2)
    library(RColorBrewer)
    library(ComplexHeatmap)
	library(sessioninfo)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_clustered.rds"))

# number of nuclei per cluster and Dataset
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Dataset)

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

pdf(file = here::here("plots", "sce_plots", "Heatmap_markers_sce.pdf"), width = 12, height = 8)

Heatmap(
  hm_mat,
  name = "z-score",
  column_title = "DG clusters mean marker expression",
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

######
# UMAP
######

# cluster labels
cluster_pops <- list(
  ExctN = c(4, 9, 1, 6, 8, 10),
  InhbN = c(19, 7, 16, 11, 20),
  Astro = c(13),
  Oligo = c(15, 12, 2, 5, 18, 14),
  Macro_Micro_T = c(17),
  OPCs = c(3),
  Endo_Mural = c(21)
    )

label_merged <- fct_collapse(colData(sce)$label,
  ExctN = as.character(cluster_pops[[1]]),
  InhbN = as.character(cluster_pops[[2]]),
  Astro = as.character(cluster_pops[[3]]),
  Oligo = as.character(cluster_pops[[4]]),
  Macro_Micro_T = as.character(cluster_pops[[5]]),
  OPCs = as.character(cluster_pops[[6]]),
  Endo_Mural = as.character(cluster_pops[[7]])
    )

label_merged <- fct_relevel(label_merged,
  c("ExctN", "InhbN", "Astro", "Oligo", "Macro_Micro_T",
      "OPCs", "Endo_Mural"))

table(label_merged)
#label_merged
#        ExctN         InhbN         Astro         Oligo Macro_Micro_T
#        60373         24090         27609         69768         13613
#         OPCs    Endo_Mural
#        12732          1090

colData(sce)$label_merged <- label_merged

colors_clusters <- list(population = c(
  ExctN = "blue",
  InhbN = "green",
  Astro = "yellow",
  Oligo = "plum3",
  Macro_Micro_T = "tan",
  OPCs = "goldenrod",
  Endo_Mural = "red3")
    )

pdf(file = here::here("plots", "sce_plots", "Merged_cluster_plot_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label_merged") +
  scale_color_manual(values = colors_clusters[[1]], name = "clusters (merged)") +
  theme_classic() +
  ggtitle("Combined HPC snRNAseq datasets clustering")

dev.off()

saveRDS(sce, file = here::here("processed-data", "sce", "sce_clustered.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
