#############################################################################
# spatial_DG_lifespan project
# Heatmap of genes associated with reactive astrocytes
# Anthony Ramnauth, June 29 2023
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

# Load gene sets

Clarke_2018 <- list(
    PAN = c("Lcn2", "Steap4", "S1pr3", "Timp1", "Hsbp1", "Cxcl10", "Cd44", "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap"),
    A1 = c("C3", "H2-T23", "Serping1", "H2-D1", "Ggta1", "Ligp1", "Gpp2", "Fbln5", "Fkbp5", "Psmb8", "Srgn", "Amigo2"),
    A2 = c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14", "Stat3")
)

# Translate from one species to the other using the orthology
PAN <- orthology[orthology$Column3 %in% Clarke_2018$PAN,]
PAN <- as.data.frame(PAN$Column1)
PAN$gene_name <- PAN$`PAN$Column1`
PAN$`PAN$Column1` <- NULL
PAN$label <- rep("PAN", 13)
A1 <- orthology[orthology$Column3 %in% Clarke_2018$A1,]
A1 <- as.data.frame(A1$Column1)
A1$gene_name <- A1$`A1$Column1`
A1$`A1$Column1` <- NULL
A1$label <- rep("A1", 22)
A2 <- orthology[orthology$Column3 %in% Clarke_2018$A2,]
A2 <- as.data.frame(A2$Column1)
A2$gene_name <- A2$`A2$Column1`
A2$`A2$Column1` <- NULL
A2$label <- rep("A2", 13)

Clarke_2018_all <- rbind(PAN, A1, A2)
Clarke_2018_all <- unique(Clarke_2018_all)

Clarke_2018_all1 <- Clarke_2018_all[! Clarke_2018_all$gene_name %in%
        setdiff(Clarke_2018_all$gene_name, rownames(spe_pseudo)),]

Clarke_heatmap <- assays(spe_pseudo)[[2]][Clarke_2018_all1$gene_name, ]
colnames(Clarke_heatmap) <- paste("logcount", 1:64, sep = "")

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

Clarke_heatmap <- scale_rows(Clarke_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Clarke_astrocytes_genemarkers_heatmap.pdf"))

Heatmap(Clarke_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(sub_type = Clarke_2018_all$label,
        col = list(sub_type = c("PAN" = "lightblue", "A1" = "dodgerblue", "A2" = "midnightblue"))),
    column_title = "Clarke et al., 2018 Gene markers for reactive astrocytes",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    row_split = Clarke_2018_all$label,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()

# Get list of gene-set from human data (Yijing Su et al., 2022) for microglia neuroinflammatory cluster
Su_astro_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Su_astro_2022.csv"))

Su_astro1_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST1",]
Su_astro6_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST6",]
Su_astro_2022 <- rbind(Su_astro1_2022, Su_astro6_2022)
Su_astro_2022 <- data.frame(gene_name = Su_astro_2022$Gene,
    label = Su_astro_2022$Cluster.ID)
Su_astro_2022 <- unique(Su_astro_2022)

Su_astro_2022 <- Su_astro_2022[! Su_astro_2022$gene_name %in%
        setdiff(Su_astro_2022$gene_name, rownames(spe_pseudo)),]

Su_heatmap <- assays(spe_pseudo)[[2]][Su_astro_2022$gene_name, ]
colnames(Su_heatmap) <- paste("logcount", 1:64, sep = "")

Su_heatmap <- scale_rows(Su_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Su_astrocytes_genemarkers_heatmap.pdf"))

Heatmap(Su_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(sub_type = Su_astro_2022$label,
        col = list(sub_type = c("AST1" = "lightblue", "AST6" = "midnightblue"))),
right_annotation = rowAnnotation(foo = anno_mark(at = c(1, 4, 73, 194, 213, 114, 212, 14, 25, 91, 112, 204, 234,
    255, 308, 669, 691, 704),
        labels = c("CD44", "GFAP", "C3", "VIM", "FKBP5", "CD109", "STAT3", "VCAN", "MALAT1",
            "MAOB", "CD38", "S100B", "AQP4", "AHDC1", "ALDH1A1", "CD81", "ILF3", "SOX9"))),
    column_title = "Su et al., 2022 Gene markers for reactive astrocytes",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    row_split = Su_astro_2022$label,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    )

dev.off()

######################################################################################################################

# Plot heatmap separately for SLM
# Relead the spe_pseudo or rename them in above code, since previously limited to DG

spe_SLM <- spe_pseudo[, which(spe_pseudo$bayesSpace_harmony_10 == "1")]
dim(spe_SLM)

## Configure column order to match age groups per BayesSpace cluster
Bayes_age_order_SLM <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6
)

## Set gene names as row names for easier plotting
rownames(spe_SLM) <- rowData(spe_SLM)$gene_name

Clarke_2018_SLM <- Clarke_2018_all[! Clarke_2018_all$gene_name %in%
        setdiff(Clarke_2018_all$gene_name, rownames(spe_SLM)),]

Clarke_SLM_heatmap <- assays(spe_SLM)[[2]][Clarke_2018_SLM$gene_name, ]
colnames(spe_SLM) <- paste("logcount", 1:16, sep = "")

Clarke_SLM_heatmap <- scale_rows(Clarke_SLM_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "SLM_Clarke_astrocytes_genemarkers_heatmap.pdf"))

Heatmap(Clarke_SLM_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = spe_SLM$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(sub_type = Clarke_2018_SLM$label,
        col = list(sub_type = c("PAN" = "lightblue", "A1" = "dodgerblue", "A2" = "midnightblue"))),
    column_title = "Clarke et al., 2018 Gene markers for reactive astrocytes",
    column_order = Bayes_age_order_SLM,
    show_column_names = FALSE,
    row_split = Clarke_2018_SLM$label,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()

Su_astro_2022_SLM <- Su_astro_2022[! Su_astro_2022$gene_name %in%
        setdiff(Su_astro_2022$gene_name, rownames(spe_SLM)),]

Su_SLM_heatmap <- assays(spe_SLM)[[2]][Su_astro_2022_SLM$gene_name, ]
colnames(spe_SLM) <- paste("logcount", 1:16, sep = "")

Su_SLM_heatmap <- scale_rows(Su_SLM_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "SLM_Su_astrocytes_genemarkers_heatmap.pdf"))

Heatmap(Su_SLM_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = spe_SLM$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(sub_type = Su_astro_2022_SLM$label,
        col = list(sub_type = c("AST1" = "lightblue", "AST6" = "midnightblue"))),
    right_annotation = rowAnnotation(foo = anno_mark(at = c(1, 4, 73, 194, 213, 114, 212, 14, 25, 91, 112, 204, 234,
    255, 308, 669, 691, 704),
        labels = c("CD44", "GFAP", "C3", "VIM", "FKBP5", "CD109", "STAT3", "VCAN", "MALAT1",
            "MAOB", "CD38", "S100B", "AQP4", "AHDC1", "ALDH1A1", "CD81", "ILF3", "SOX9"))),
    column_title = "Su et al., 2022 Gene markers for reactive astrocytes",
    column_order = Bayes_age_order_SLM,
    show_column_names = FALSE,
    row_split = Su_astro_2022_SLM$label,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    )

dev.off()

