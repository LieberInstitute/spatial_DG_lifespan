#####################################################
# spatial_DG_lifespan project
# Sub-clustering of sce object for cells
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
    library(forcats)
    library(RColorBrewer)
    library(dplyr)
    library(tidyr)
    library(harmony)
    library(ComplexHeatmap)
    library(sessioninfo)
})

# load saved sce object

sce <-
    readRDS(here::here("processed-data", "sce", "sce_clustered.rds"))

# Add character to sub-cluster columns to denote cell type for each cluster number

colData(sce)$label_excitatory <-
    sub("^", "ExctN", colData(sce)$label_excitatory)
colData(sce)$label_inhibitory <-
    sub("^", "IntN", colData(sce)$label_inhibitory)
colData(sce)$label_astro <-
    sub("^", "Astro", colData(sce)$label_astro)
colData(sce)$label_oligo <-
    sub("^", "Oligo", colData(sce)$label_oligo)
colData(sce)$label_immune <-
    sub("^", "Immun", colData(sce)$label_immune)
colData(sce)$label_OPCs <- sub("^", "OPC", colData(sce)$label_OPCs)
colData(sce)$label_COP <- sub("^", "COP", colData(sce)$label_COP)
colData(sce)$label_Endo_Mural <-
    sub("^", "Endo_Mur", colData(sce)$label_Endo_Mural)

# Combine columns into one column for cell-type subclusters

dfc <- data.frame(
    label_excitatory = colData(sce)$label_excitatory,
    label_inhibitory = colData(sce)$label_inhibitory,
    label_astro = colData(sce)$label_astro,
    label_oligo = colData(sce)$label_oligo,
    label_immune = colData(sce)$label_immune,
    label_OPCs = colData(sce)$label_OPCs,
    label_COP = colData(sce)$label_COP,
    label_Endo_Mural = colData(sce)$label_Endo_Mural
)

dfc <- dfc %>%
    mutate(
        cell_type = coalesce(
            label_excitatory,
            label_inhibitory,
            label_astro,
            label_oligo,
            label_immune,
            label_OPCs,
            label_COP,
            label_Endo_Mural
        )
    )

colData(sce)$cell_type <- dfc$cell_type

# Markers chosen from the publications of each dataset

markers <- c(
    ## Neurons
    "RBFOX3",
    "SNAP25",
    "SYT1",
    ## Excitatory Neurons
    "SLC17A7",
    # GC
    "PROX1",
    # im GC
    "DCX",
    "BHLHE22",
    "STMN1",
    # Mossy Cells
    "ARHGAP24",
    "DLC1",
    # CA3 PNs
    "CFAP299",
    "SYN3",
    # CA2 PNs
    "HGF",
    # CA1 PNs
    "ACVR1C",
    "SYT13",
    # Sub PNs
    "ROBO1",
    "COL5A2",
    # Cajal–Retzius
    "RELN",
    ## Inhibitory Neurons
    "GAD1",
    "GAD2",
    # inhibitory subpopulations (Some from Lukas LC & Keri Martinowich 2022-07-22)
    "SST",
    "KIT",
    "CALB1",
    "CALB2",
    "TAC1",
    "CNR1",
    "PVALB",
    "CORT",
    "VIP",
    "NPY",
    "CRHBP",
    "CCK",
    "HTR3A",
    "NR2F2",
    "LAMP5",
    # Astrocytes
    "AQP4",
    "GFAP",
    "CHRDL1",
    # Oligodendrocytes
    "MOBP",
    # macrophages / microglia
    "CD163",
    "C3",
    "PTPRC",
    "C1QB",
    # OPCs
    "PDGFRA",
    "VCAN",
    # COP
    "GPR17",
    "ADAM33",
    # endothelial / mural (RBPMS)
    "CLDN5",
    "FLT1",
    "RBPMS",
    # T cells
    "SKAP1",
    "CD247",
    # Progenitors
    "PAX6",
    "HOPX",
    "EOMES"
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
colors_markers <- list(
    marker = c(
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
    )
)

marker_labels <-
    factor(marker_labels, levels = unique(marker_labels))

# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x)
    split(seq(along = x), x)

cell_idx <- splitit(colData(sce)$cell_type)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$SYMBOL

hm_mat <-
    t(do.call(cbind, lapply(cell_idx, function(i)
        rowMeans(dat[markers, i]))))

# convert to z-scores
scale_rows = function(x) {
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

# set the row order for subclusters
roword <- c("ExctN1", "ExctN2", "ExctN3", "ExctN4", "ExctN5", "ExctN6", "ExctN7",
    "ExctN8", "ExctN9",
    "IntN1", "IntN2", "IntN3", "IntN4", "IntN5", "IntN6", "IntN7",
    "Astro1", "Astro2", "Astro3", "Astro4", "Astro5", "Astro6",
    "Oligo1", "Oligo2", "Oligo3", "Oligo4",
    "Immun1", "Immun2", "Immun3", "Immun4", "Immun5",
    "OPC1",
    "COP1",
    "Endo_Mur1", "Endo_Mur2")

pdf(
    file = here::here("plots", "sce_plots", "Heatmap_markers_sce.pdf"),
    width = 12,
    height = 8
)

Heatmap(
    hm_mat,
    name = "z-score",
    column_title = "DG clusters markers",
    column_title_gp = gpar(fontface = "bold"),
    bottom_annotation = col_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = roword,
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
  ExctN = c("ExctN1", "ExctN2", "ExctN3", "ExctN4", "ExctN5", "ExctN6", "ExctN7",
    "ExctN8", "ExctN9"),
  IntN = c("IntN1", "IntN2", "IntN3", "IntN4", "IntN5", "IntN6", "IntN7"),
  Astro = c("Astro1", "Astro2", "Astro3", "Astro4", "Astro5", "Astro6"),
  Oligo = c("Oligo1", "Oligo2", "Oligo3", "Oligo4"),
  Immune = c("Immun1", "Immun2", "Immun3", "Immun4", "Immun5"),
  OPCs = c("OPC1"),
  COP = c("COP1"),
  Endo_Mural = c("Endo_Mur1", "Endo_Mur2")
    )

label_merged <- fct_collapse(colData(sce)$cell_type,
  ExctN = as.character(cluster_pops[[1]]),
  IntN = as.character(cluster_pops[[2]]),
  Astro = as.character(cluster_pops[[3]]),
  Oligo = as.character(cluster_pops[[4]]),
  Immune = as.character(cluster_pops[[5]]),
  OPCs = as.character(cluster_pops[[6]]),
  COP = as.character(cluster_pops[[7]]),
  Endo_Mural = as.character(cluster_pops[[8]])
    )

label_merged <- fct_relevel(label_merged,
  c("ExctN", "IntN", "Astro", "Oligo", "Immune",
      "OPCs", "COP", "Endo_Mural"))

table(label_merged)
#label_merged
#     ExctN      IntN      Astro      Oligo     Immune       OPCs        COP
#     60407      22288      27609      69768      13613      12172        560
#Endo_Mural
#      2858

colData(sce)$label_merged <- label_merged

colors_clusters <- list(population = c(
  ExctN = "blue",
  IntN = "green",
  Astro = "yellow",
  Oligo = "plum3",
  Immune = "tan",
  OPCs = "goldenrod",
  COP = "khaki",
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
