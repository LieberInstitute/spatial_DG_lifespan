###############################
# spatial_DG_lifespan project
# PCA on neuroblast 2 gene set
# Anthony Ramnauth, Apr 21 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(ggplot2)
    library(ggnewscale)
    library(PCAtools)
    library(spatialLIBD)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

## Set gene names as row names for easier orthology conversion
rownames(spe) <- rowData(spe)$gene_name

# order spe observations according to age
spe <- spe[, order(spe$age)]

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

setdiff(Hochgerner_2018_NB2$gene_names, rownames(spe))

Hochgerner_2018_NB2 <- Hochgerner_2018_NB2[! Hochgerner_2018_NB2$gene_names %in%
        setdiff(Hochgerner_2018_NB2$gene_names, rownames(spe)),]

# Dimensionality reduction on NB2 geneset
set.seed(12345)
spe <- runPCA(spe, subset_row = Hochgerner_2018_NB2$gene_names, ncomponents = 10)
set.seed(12345)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")


# Plot UMAP before harmony
pdf(file = here::here("plots", "gene_set_enrichment", "NB2_UMAP.pdf"))

ggplot(data.frame(reducedDim(spe, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = spe$bayesSpace_harmony_10,
        alpha = 0.5
    )) +
    geom_point() +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace") +
    theme_bw() +
    theme_classic()

dev.off()

pdf(file = here::here("plots", "gene_set_enrichment", "NB2_PC1vsPC2.pdf"))

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "bayesSpace_harmony_10",
    point_size = 2,
    point_alpha = 0.5) +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace") +
    theme_bw() +
    theme_classic()

plotReducedDim(
    spe,
    dimred = "PCA",
    ncomponents = 2,
    colour_by = "age",
    point_size = 2,
    point_alpha = 0.5) +
    labs(color = "age") +
    theme_bw() +
    theme_classic()

dev.off()

spe$PCA1 <- reducedDims(spe)$PCA[,1]
spe$PCA2 <- -(reducedDims(spe)$PCA[,2])

vis_grid_gene(spe = spe,
    geneid = "PCA1",
    spatial = TRUE,
    viridis - FALSE,
    pdf = here::here("plots", "gene_set_enrichment", "NB2_PC1.pdf"),
    image_id = "lowres",
    alpha = 0.5,
    point_size = 2
)

vis_grid_gene(spe = spe,
    geneid = "PCA2",
    spatial = TRUE,
    viridis - FALSE,
    pdf = here::here("plots", "gene_set_enrichment", "NB2_PC2.pdf"),
    image_id = "lowres",
    alpha = 0.5,
    point_size = 2
)
