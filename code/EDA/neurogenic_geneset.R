#############################################################################
# spatial_DG_lifespan project
# Heatmap of genes associated with neurogenic, & progenitor cells
# Anthony Ramnauth, June 05 2023
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
    library(sessioninfo)
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

# Translate from one species to the other using the orthology
Hochgerner_2018_nIPC <- orthology[orthology$Column3 %in% Hochgerner_2018$nIPC,]
Hochgerner_2018_nIPC <- Hochgerner_2018_nIPC$Column1
Hochgerner_2018_nIPC <- as.data.frame(Hochgerner_2018_nIPC)
Hochgerner_2018_nIPC$label <- rep("nIPC", 113)
names(Hochgerner_2018_nIPC)[names(Hochgerner_2018_nIPC) == "Hochgerner_2018_nIPC"] <- "gene_names"

Hochgerner_2018_NB1 <- orthology[orthology$Column3 %in% Hochgerner_2018$Neuroblast1,]
Hochgerner_2018_NB1 <- Hochgerner_2018_NB1$Column1
Hochgerner_2018_NB1 <- as.data.frame(Hochgerner_2018_NB1)
Hochgerner_2018_NB1$label <- rep("NB1", 41)
names(Hochgerner_2018_NB1)[names(Hochgerner_2018_NB1) == "Hochgerner_2018_NB1"] <- "gene_names"

Hochgerner_2018_NB2 <- orthology[orthology$Column3 %in% Hochgerner_2018$Neuroblast2,]
Hochgerner_2018_NB2 <- Hochgerner_2018_NB2$Column1
Hochgerner_2018_NB2 <- as.data.frame(Hochgerner_2018_NB2)
Hochgerner_2018_NB2$label <- rep("NB2", 92)
names(Hochgerner_2018_NB2)[names(Hochgerner_2018_NB2) == "Hochgerner_2018_NB2"] <- "gene_names"

Hochgerner_2018_all <- rbind(Hochgerner_2018_nIPC, Hochgerner_2018_NB1, Hochgerner_2018_NB2)

Hochgerner_2018_all <- Hochgerner_2018_all[! Hochgerner_2018_all$gene_names %in%
        setdiff(Hochgerner_2018_all$gene_names, rownames(spe_pseudo)),]

Hochgerner_heatmap <- assays(spe_pseudo)[[2]][Hochgerner_2018_all$gene_names, ]
colnames(Hochgerner_heatmap) <- paste("logcount", 1:64, sep = "")

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

Hochgerner_heatmap <- scale_rows(Hochgerner_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "nIPCs_NBs_genemarkers_heatmap.pdf"),
    width = 12, height = 14)

Heatmap(Hochgerner_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(cell_marker = Hochgerner_2018_all$label),
    column_title = "Gene markers for nIPCs & NBs",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    row_order = Hochgerner_2018_all$gene_names
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

conserved_heatmap <- assays(spe_pseudo)[[2]][conserved_final, ]
colnames(conserved_heatmap) <- paste("logcount", 1:64, sep = "")

conserved_heatmap <- scale_rows(conserved_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "mouse_macaque_conserved_genemarkers_heatmap.pdf"),
    width = 12, height = 14)

Heatmap(conserved_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Neurogenic markers conserved between mouse and macaque",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    )

dev.off()

# Get list of imGC gene-set from human-lifespan data (Yi Zhou et al., 2022)
Zhou_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Zhou_2022.csv"))
Zhou_2022$Number <- NULL

common_imGC <- c(Zhou_2022$Common.genes[1:84], Zhou_2022$Human.specific.genes[1:76])

common_imGC_final <- common_imGC[! common_imGC %in%
        setdiff(common_imGC, rownames(spe_pseudo))]

common_imGC_heatmap <- assays(spe_pseudo)[[2]][common_imGC_final, ]
colnames(common_imGC_heatmap) <- paste("logcount", 1:64, sep = "")

common_imGC_heatmap <- scale_rows(common_imGC_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Common_imGC_genemarkers_heatmap.pdf"),
    width = 12, height = 14)

Heatmap(common_imGC_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "imGC markers conserved between mouse and human",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    )

dev.off()

