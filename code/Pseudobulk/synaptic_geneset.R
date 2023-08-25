#############################################################################
# spatial_DG_lifespan project
# Heatmap of genes associated with neuropil
# Anthony Ramnauth, July 04 2023
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
    library(circlize)
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

# Use table from Jax labs for mouse human orthology (gene names are same in rats)
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Get list of gene-set from mouse data (Iván J. Cajigas et al., 2012) for synaptically enriched gene sets
Cajigas_2012 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Cajigas_2012.csv"))

# Translate from one species to the other using the orthology
Cajigas_2012_rat <- orthology[orthology$Column3 %in% Cajigas_2012$Gene.Symbol,]
Cajigas_2012_rat <- as.data.frame(Cajigas_2012_rat$Column1)
colnames(Cajigas_2012_rat)[1] <- "gene_name"

# Get list of gene-set from mouse data (Robert R. Stickels et al., 2020) for dendritically enriched gene sets
Stickels_2020 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Stickels_2020.csv"))

# Translate from one species to the other using the orthology
Dcluster_1 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.1,]
Dcluster_1 <- as.data.frame(Dcluster_1$Column1)
Dcluster_1$label <- rep("cluster_1", 62)
colnames(Dcluster_1)[1] <- "gene_name"

Dcluster_2 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.2,]
Dcluster_2 <- as.data.frame(Dcluster_2$Column1)
Dcluster_2$label <- rep("cluster_2", 55)
colnames(Dcluster_2)[1] <- "gene_name"

Dcluster_3 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster.3,]
Dcluster_3 <- as.data.frame(Dcluster_3$Column1)
Dcluster_3$label <- rep("cluster_3", 34)
colnames(Dcluster_3)[1] <- "gene_name"

Dcluster_4 <- orthology[orthology$Column3 %in% Stickels_2020$Cluster4,]
Dcluster_4 <- as.data.frame(Dcluster_4$Column1)
Dcluster_4$label <- rep("cluster_4", 57)
colnames(Dcluster_4)[1] <- "gene_name"

Dcluster_all <- rbind(Dcluster_1, Dcluster_2, Dcluster_3, Dcluster_4)

# Get list of gene-set from mouse data (Muchun Niu et al., 2023) for synaptically enriched gene sets
Niu_2023 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Niu_2023.csv"))

# Subset for Synapse_ExDG and Synapse_In
Niu_2023_DG <- Niu_2023[Niu_2023$cluster == "Synapse_ExDG" | Niu_2023$cluster == "Synapse_In",]

#######################################################################################################

# Find each geneset in spe_pseudo data and plot

Cajigas_2012_rat <- Cajigas_2012_rat[! Cajigas_2012_rat$gene_name %in%
        setdiff(Cajigas_2012_rat$gene_name, rownames(spe_pseudo)),]

Cajigas_heatmap <- assays(spe_pseudo)[[2]][Cajigas_2012_rat, ]
colnames(Cajigas_heatmap) <- paste("logcount", 1:64, sep = "")

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

Cajigas_heatmap <- scale_rows(Cajigas_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Rat_synaptome_heatmap.pdf"))

Heatmap(Cajigas_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = "Cajigas et al., 2012 Gene markers for rat synaptome",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    )

dev.off()

#########################################################################################################

# Find each geneset in spe_pseudo data and plot

Dcluster <- Dcluster_all[! Dcluster_all$gene_name %in%
        setdiff(Dcluster_all$gene_name, rownames(spe_pseudo)),]

Stickels_heatmap <- assays(spe_pseudo)[[2]][Dcluster$gene_name, ]
colnames(Stickels_heatmap) <- paste("logcount", 1:64, sep = "")

Stickels_heatmap <- scale_rows(Stickels_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Mouse_dendritic_heatmap.pdf"))

Heatmap(Stickels_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(cluster = Dcluster$label,
        col = list(cluster = c("cluster_1" = "lightblue", "cluster_2" = "dodgerblue", "cluster_3" = "midnightblue", "cluster_4" = "black"))),
    column_title = "Stickels et al., 2020 mouse dendritically enriched gene markers",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    row_split = Dcluster$label,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    )

dev.off()

###############################################################################################################


# Find each geneset in spe_pseudo data and plot

Niu_2023_DG <- Niu_2023_DG[! Niu_2023_DG$gene %in%
        setdiff(Niu_2023_DG$gene, rownames(spe_pseudo)),]

Niu_heatmap <- assays(spe_pseudo)[[2]][Niu_2023_DG$gene, ]
colnames(Niu_heatmap) <- paste("logcount", 1:64, sep = "")

Niu_heatmap <- scale_rows(Niu_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Human_synaptome_heatmap.pdf"))

Heatmap(Niu_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(cluster = Niu_2023_DG$cluster,
        col = list(cluster = c("Synapse_ExDG" = "lightblue", "Synapse_In" = "midnightblue"))),
    column_title = "Niu et al., 2023 human DG synaptome gene markers",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    row_split = Niu_2023_DG$cluster,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    )

dev.off()

#################################################################################################################

# Isolate ML enriched gene sets and plot only in ML

Niu_2023_ExDG <- Niu_2023[Niu_2023$cluster == "Synapse_ExDG",]

neuropil_ML <- unique(c(Dcluster_3$gene_name, Dcluster_4$gene_name, Niu_2023_ExDG$gene))

# Limit pseudo-bulk data to Molecular layer
pseudo_ML <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "2")]
dim(pseudo_ML)

# Find genes only in dataset
neuropil_ML <- neuropil_ML[! neuropil_ML %in%
        setdiff(neuropil_ML, rownames(spe_pseudo))]

## Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6
)

# Add logcounts for all neuropil genes
neuropil_heatmap <- assays(pseudo_ML)[[2]][neuropil_ML, ]
colnames(neuropil_heatmap) <- paste("logcount", 1:16, sep = "")


# Plot heatmap
pdf(file = here::here("plots", "pseudobulked", "ML_enriched_synaptome_heatmap.pdf"))

col_fun = colorRamp2(c(0, 5, 15), c("blue", "white", "red"))

Heatmap(neuropil_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = pseudo_ML$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    right_annotation = rowAnnotation(foo = anno_mark(at = c(171, 137, 203, 264, 122, 83, 186, 156, 223, 229, 149, 78,
        24, 25, 50, 73, 77, 27, 80, 74, 76, 72, 2, 220, 16, 261, 231, 112),
        labels = c("SNAP25", "SORT1", "CNKSR2", "ATP6V1B2", "PEA15", "SOWAHA", "MAST3", "PSAP", "BMERB1", "AK5", "APP", "RPS3",
            "RPL13", "RPL13A", "EEF1G", "RPLP1", "RPS24", "RPS15", "RPS5", "RPS11", "RPS2", "RPLP0", "ARHGAP5", "NRXN1", "LHX2",
            "PRRC2C", "BCL11B", "RFX3"))),
    col = col_fun,
    column_title = "Synaptic/dendrtic enriched gene markers in ML",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = TRUE
    )

dev.off()
