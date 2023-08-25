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
    library(circlize)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on BayesSpace clusters for DG
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("1", "2", "4", "6", "7")]

bayes_df <- data.frame(spe_pseudo$BayesSpace)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("1", spe_pseudo.BayesSpace) ~ "SLM",
        grepl("2", spe_pseudo.BayesSpace) ~ "ML",
        grepl("4", spe_pseudo.BayesSpace) ~ "CA3&4",
        grepl("6", spe_pseudo.BayesSpace) ~ "SGZ",
        grepl("7", spe_pseudo.BayesSpace) ~ "GCL"
    ))

colData(spe_pseudo)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("SLM", "ML", "CA3&4", "SGZ", "GCL"))

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 27, 18, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 43, 34, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 59, 50, 57, 58, 62, 52, 51, 53, 55, 54,
    80, 76, 77, 79, 72, 65, 75, 66, 73, 74, 78, 68, 67, 69, 71, 70
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
colnames(Clarke_heatmap) <- paste("logcount", 1:80, sep = "")

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Clarke_astrocytes_genemarkers_heatmap.pdf"))

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

Heatmap(Clarke_heatmap,
    name = "mean\nnorm logcounts",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(spatial_domain = c("SLM" = "black", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    left_annotation = rowAnnotation(sub_type = Clarke_2018_all1$label,
        col = list(sub_type = c("PAN" = "lightblue", "A1" = "dodgerblue", "A2" = "midnightblue"))),
    col = col_fun,
    column_title = "Clarke et al., 2018 Gene markers for reactive astrocytes",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
    row_split = Clarke_2018_all1$label,
    show_row_names = TRUE,
    cluster_rows = TRUE,
    )

dev.off()

# Get list of gene-set from human data (Yijing Su et al., 2022) for microglia neuroinflammatory cluster
Su_astro_2022 <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Su_astro_2022.csv"))

Su_astro6_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST6",]

Su_astro_2022 <- Su_astro6_2022[! Su_astro6_2022$Gene %in%
        setdiff(Su_astro6_2022$Gene, rownames(spe_pseudo)),]

Su_heatmap <- assays(spe_pseudo)[[2]][Su_astro_2022$Gene, ]
colnames(Su_heatmap) <- paste("logcount", 1:64, sep = "")

Su_heatmap <- scale_rows(Su_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "Su_astrocytes_genemarkers_heatmap.pdf"),
    width = 7, height = 9)

Heatmap(Su_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(BayesSpace = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(BayesSpace = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
right_annotation = rowAnnotation(foo = anno_mark(at = c(213, 362, 225, 143, 369, 43, 44,
    22, 101, 40, 345, 289, 52, 395, 280, 201, 17, 103, 293, 231, 37, 331, 25, 157, 185,
    245, 142, 279, 68, 89, 164, 121, 5, 155, 145, 167, 48, 141, 81, 371, 15, 304, 107,
    133, 188, 83),
        labels = c("C6orf89", "ICE1", "ZNF106", "FAIM2", "TRIM37", "PPP3R1", "PEA15", "IDS", "ATP1A1",
            "NCDN", "SMIM14", "CDS2", "PEBP1", "VAPA", "JAK1", "ENO2", "SNAP25", "NAPB", "ATP6V1C1",
            "ATP6V1E1", "VSNL1", "AGPAT3", "PACSIN1", "ZNF483", "TMX4", "CYFIP2", "PSAP", "FRMPD4",
            "TPPP", "FRRS1L", "TGOLN2", "TMEM30A", "CREG2", "PREPL", "PJA2", "CLSTN1", "TOMM20",
            "SERINC3", "DNM1", "ATP6V0A1", "SLC17A7", "CD99L2", "DYNC1H1", "CHN1", "SYNJ1", "MTURN"))),
    column_title = "Su et al., 2022 Gene markers for AST6 cluster",
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    column_split = spe_pseudo$BayesSpace,
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
colnames(Clarke_SLM_heatmap) <- paste("logcount", 1:16, sep = "")

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

Su_astro1_2022 <- Su_astro_2022[Su_astro_2022$Cluster.ID == "AST1",]

Su_astro_2022_SLM <- Su_astro1_2022[! Su_astro1_2022$Gene %in%
        setdiff(Su_astro1_2022$Gene, rownames(spe_SLM)),]

Su_SLM_heatmap <- assays(spe_SLM)[[2]][Su_astro_2022_SLM$Gene, ]
colnames(Su_SLM_heatmap) <- paste("logcount", 1:16, sep = "")

Su_SLM_heatmap <- scale_rows(Su_SLM_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked", "SLM_Su_astrocytes_genemarkers_heatmap.pdf"))

Heatmap(Su_SLM_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(age = spe_SLM$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    right_annotation = rowAnnotation(foo = anno_mark(at = c(49, 44, 91, 187, 4, 34, 242,
        28, 149, 197, 328, 330, 158, 82, 204, 33),
        labels = c("ITGB4", "ITPKB", "MAOB", "SOD2", "GFAP", "SPARC", "WDR11", "SNED1", "EFEMP1", "SLC9A9",
            "GSN", "TRIM2", "AGFG2", "PFKFB2", "S100B", "USH1C"))),
    column_title = "Su et al., 2022 Gene markers for AST1 cluster",
    column_order = Bayes_age_order_SLM,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    )

dev.off()

