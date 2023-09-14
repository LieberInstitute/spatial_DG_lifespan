###########################################
# spatial_DG_lifespan project
# Reduced Dimensions for Franjic et al 2022
# Anthony Ramnauth, June 13 2022
###########################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(ggplot2)
    library(PCAtools)
    library(ggnewscale)
    library(spatialLIBD)
    library(schex)
})

# Load SCE
sce <- readRDS(file = here::here("processed-data", "sce", "sce_sestan_DG_final.rds"))

dim(sce)

# Feature selection
dec <- modelGeneVar(sce, block = sce$sample_id)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
head(top_hvgs)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
head(top_hvgs)

# Dimensionality reduction
set.seed(12345)
sce <- runPCA(sce, subset_row = top_hvgs, ncomponents = 50)
set.seed(12345)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- c("UMAP1", "UMAP2")

cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3",
    "VLMC" = "firebrick1", "Pericyte" = "red2", "Endoth" = "red", "SMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "yellow", "Astro_2" = "yellow3")

# Plot UMAP before harmony
pdf(file = here::here("plots", "Cell_Type_Deconvolution", "sestan_UMAP_sce.pdf"))
ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce$Cell_Type),
        alpha = 0.5
    )) +
    scale_colour_manual(values = cell_colors) +
    geom_point() +
    labs(color = "Cell Type") +
    theme_bw() +
    theme_classic()

dev.off()
