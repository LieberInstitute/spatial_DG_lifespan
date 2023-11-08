#############################################################################
# spatial_DG_lifespan project
# Heatmap of genes associated NB2
# Anthony Ramnauth, Aug 22 2023
#############################################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on BayesSpace clusters for DG
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("2", "4", "6", "7")]

bayes_df <- data.frame(spe_pseudo$BayesSpace)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("2", spe_pseudo.BayesSpace) ~ "ML",
        grepl("4", spe_pseudo.BayesSpace) ~ "CA3&4",
        grepl("6", spe_pseudo.BayesSpace) ~ "SGZ",
        grepl("7", spe_pseudo.BayesSpace) ~ "GCL"
    ))

colData(spe_pseudo)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("ML", "CA3&4", "SGZ", "GCL"))

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 27, 18, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 43, 34, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 59, 50, 57, 58, 62, 52, 51, 53, 55, 54
)

## Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Get list of neurogenesis gene-set from mouse data (Hochgerner et al., 2018)
Hochgerner_2018 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Hochgerner_2018.csv"))

Hochgerner_2018_NB2 <- orthology[orthology$Column3 %in% Hochgerner_2018$Neuroblast2,]
Hochgerner_2018_NB2 <- Hochgerner_2018_NB2$Column1
Hochgerner_2018_NB2 <- as.data.frame(Hochgerner_2018_NB2)
Hochgerner_2018_NB2$label <- rep("NB2", 92)
names(Hochgerner_2018_NB2)[names(Hochgerner_2018_NB2) == "Hochgerner_2018_NB2"] <- "gene_names"

Hochgerner_2018_NB2$gene_names <- unique(Hochgerner_2018_NB2$gene_names)

setdiff(Hochgerner_2018_NB2$gene_names, rownames(spe_pseudo))

Hochgerner_2018_NB2 <- Hochgerner_2018_NB2[! Hochgerner_2018_NB2$gene_names %in%
        setdiff(Hochgerner_2018_NB2$gene_names, rownames(spe_pseudo)),]

Hochgerner_heatmap <- assays(spe_pseudo)[[2]][Hochgerner_2018_NB2$gene_names, ]
colnames(Hochgerner_heatmap) <- paste("logcount", 1:64, sep = "")

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "NB2_genemarkers_heatmap.pdf"), width = 7, height = 12.5)

Heatmap(Hochgerner_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Hochgerner et al., 2018 NB2 markers",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

Hochgerner_heatmap <- scale_rows(Hochgerner_heatmap)

pdf(file = here::here("plots", "pseudobulked", "zscores_NB2_genemarkers_heatmap.pdf"), width = 7, height = 12.5)

Heatmap(Hochgerner_heatmap,
    name = "mean\nlog norm counts\n(centered & scaled)",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Hochgerner et al., 2018 NB2 markers",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()

# List taken from Hao et al 2022 for monkey and macaque
conserved_macaque_mice_neurogenic <- c(
    "ID4", "TMEM47", "MLC1", "ALDOC", "SLC1A3", "SLC1A2", "SOX2", "MKI67", "CKS2",
    "CENPF", "SMC4", "TOP2A", "TMPO", "FXYD6", "NNAT", "SOX4", "DPYSL3", "STMN2",
    "CALB2", "SEMA3C"
)

conserved_final <- conserved_macaque_mice_neurogenic[! conserved_macaque_mice_neurogenic %in%
        setdiff(conserved_macaque_mice_neurogenic, rownames(spe_pseudo))]

conserved_final <- as.data.frame(conserved_final)
conserved_final$gene_name <- conserved_final$conserved_final
conserved_final$conserved_final <- NULL
conserved_final$label <- rep("conserved", 17)

unique_macaque <- c(
    "FGFR3", "PON2", "F3", "ETNPPL", "WIF1", "ZFP36", "ATP14A", "RGCC", "TIMP3", "PMP22", "SERPINE2", "FNBP1L",
    "DPYSL5", "IGFBP2", "CD9", "SLC35F1", "CCND1", "CACNG4", "ASCL1", "SULF2", "ANXA5", "VCAN"
)

unique_macaque_final <- unique_macaque[! unique_macaque %in%
        setdiff(unique_macaque, rownames(spe_pseudo))]

unique_macaque_final <- as.data.frame(unique_macaque_final)
unique_macaque_final$gene_name <- unique_macaque_final$unique_macaque_final
unique_macaque_final$unique_macaque_final <- NULL
unique_macaque_final$label <- rep("macaque.only", 21)

total_mouse_macaque <- rbind(conserved_final, unique_macaque_final)

conserved_heatmap <- assays(spe_pseudo)[[2]][total_mouse_macaque$gene_name, ]
colnames(conserved_heatmap) <- paste("logcount", 1:64, sep = "")

conserved_heatmap <- scale_rows(conserved_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "mouse_macaque_conserved_genemarkers_heatmap.pdf"),
    width = 8, height = 6)

Heatmap(conserved_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(species = total_mouse_macaque$label,
        col = list(species = c("conserved" = "lightblue", "macaque.only" = "midnightblue"))),
    column_title = "Hao et al., 2022 neurogenic markers between mouse and macaque",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    row_split = total_mouse_macaque$label,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()

# Get list of imGC gene-set from human-lifespan data (Yi Zhou et al., 2022)
Zhou_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Zhou_2022.csv"))
Zhou_2022$Number <- NULL

common_imGC <- data.frame(Zhou_2022$Common.genes[1:84],
    gene_names = Zhou_2022$Common.genes[1:84])
common_imGC$label <- rep("Common", 84)
common_imGC$Zhou_2022.Common.genes.1.84. <- NULL

human_specific_imGC <- data.frame(Zhou_2022$Human.specific.genes[1:76],
    gene_names = Zhou_2022$Human.specific.genes[1:76])
human_specific_imGC$label <- rep("Human.Specific", 76)
human_specific_imGC$Zhou_2022.Human.specific.genes.1.76. <- NULL

total_imGC <- rbind(common_imGC, human_specific_imGC)

total_imGC_final <- total_imGC[! total_imGC$gene_names %in%
        setdiff(total_imGC$gene_names, rownames(spe_pseudo)),]

total_imGC_heatmap <- assays(spe_pseudo)[[2]][total_imGC_final$gene_names, ]
colnames(total_imGC_heatmap) <- paste("logcount", 1:64, sep = "")

total_imGC_heatmap <- scale_rows(total_imGC_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Common_imGC_genemarkers_heatmap.pdf"),
    width = 8.5, height = 9)

Heatmap(total_imGC_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(species = total_imGC_final$label,
        col = list(species = c("Common" = "lightblue", "Human.Specific" = "midnightblue"))),
    right_annotation =
        rowAnnotation(foo = anno_mark(at = c(1, 3, 17, 21, 24, 25, 43, 48, 49, 50, 57, 61, 72, 73, 79, 115, 125, 136, 142),
        labels = c("ATP11C", "BHLHE22", "DPYSL5", "FGF13", "FXYD6", "FXYD7", "MARCKSL1", "NEUROD1", "NEUROD2", "NEUROD6", "PPP1R14C",
            "PROX1", "STMN3", "STMN4", "TUBA1A", "KCNQ5", "PCDH20", "ROBO2", "SOX5"))),
    column_title = "Zhou et al., 2022 imGC markers shared with mouse or unique in human",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    row_split = total_imGC_final$label
    )

dev.off()

#############################################################################################################
# Combine all gene markers for one big heatmap

all_neurogenic <- unique(c(Hochgerner_2018_all$gene_names, total_mouse_macaque$gene_name, total_imGC_final$gene_names))


all_neurogenic_heatmap <- assays(spe_pseudo)[[2]][all_neurogenic, ]
colnames(all_neurogenic_heatmap) <- paste("logcount", 1:64, sep = "")

all_neurogenic_heatmap <- scale_rows(all_neurogenic_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Combined_neurogenic_heatmap.pdf"),
    width = 8.5, height = 9)

Heatmap(all_neurogenic_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    right_annotation =
        rowAnnotation(foo = anno_mark(at = c(101, 346, 303, 139, 203, 282, 126, 189, 185, 242,
            30, 238, 201, 195, 57, 138, 23, 89, 243, 323, 15, 233, 245, 322, 90, 326, 263, 261,
            266, 222, 170, 314, 299, 325, 12, 275, 210, 329, 140, 254, 224, 274, 273, 316, 225, 330, 309),
        labels = c("ARHGAP15", "WIPF3", "HPCA", "KLHDC8A", "IGFBP2", "ASIC2", "FRMD3", "SOX4", "SMC4",
            "LRRN3", "NEUROG2", "IVNS1ABP", "FNBP1L", "WIF1", "VCAN", "KIF26B", "LIMD1", "SOX11", "MARCKSL1",
            "PDZRN3", "DLL3", "HEG1", "NCAN", "PDLIM5", "SOX5", "PTPRD", "SCRN1", "RPS2", "SMARCD3",
            "EEF2", "SSTR2", "MGAT5", "GFRA1", "PMEPA1", "CKS2", "TUBB2A", "BZW1", "RASGRF1",
            "LIN7A", "PPP1R14C", "FGF13", "TUBA1A", "TSPAN18", "NCALD", "FGFR1", "RGS8", "KCNAB1"))),
    column_title = "Combined neurogenic markers",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = FALSE,
    cluster_rows = TRUE
    )

dev.off()
