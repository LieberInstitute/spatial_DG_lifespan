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

pdf(
    file = here::here("plots", "sce_plots", "Heatmap_markers_sce.pdf"),
    width = 12,
    height = 8
)

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

colors_clusters <- list(
    population = c(
        ExctN1 = "dodgerblue",
        ExctN2 = "blue",
        ExctN3 = "dodgerblue1",
        ExctN4 = "blue1",
        ExctN5 = "dodgerblue2",
        ExctN6 = "blue2",
        ExctN7 = "dodgerblue3",
        ExctN8 = "blue3",
        ExctN9 = "dodgerblue4",
        ExctN10 = "blue4",
        IntN1 = "green",
        IntN2 = "seagreen",
        IntN3 = "green1",
        IntN4 = "seagreen1",
        IntN5 = "green2",
        IntN6 = "seagreen2",
        IntN7 = "green3",
        IntN8 = "seagreen3",
        IntN9 = "green4",
        Astro1 = "yellow",
        Astro2 = "gold",
        Astro3 = "yellow1",
        Astro4 = "gold1",
        Astro5 = "yellow2",
        Astro6 = "gold2",
        Oligo1 = "plum",
        Oligo2 = "plum1",
        Oligo3 = "plum2",
        Oligo4 = "plum3",
        Immun1 = "tan",
        Immun2 = "tan1",
        Immun3 = "tan2",
        Immun4 = "tan3",
        Immun5 = "tan4",
        OPC1 = "purple",
        OPC2 = "purple1",
        OPC3 = "purple2",
        OPC4 = "purple3",
        OPC5 = "purple4",
        Endo_Mur1 = "red",
        Endo_Mur2 = "red1",
        Endo_Mur3 = "red2"
    )
)

pdf(file = here::here("plots", "sce_plots", "Temp_celltype_cluster_plot_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "cell_type") +
    scale_color_manual(values = colors_clusters[[1]], name = "clusters") +
    theme_classic() +
    ggtitle("Combined HPC snRNAseq datasets intermediate cell-type clustering")

dev.off()

# Save new sce object
saveRDS(sce,
    file = here::here("processed-data", "sce", "sce_clustered.rds"))
