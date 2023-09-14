###############################
# spatial_DG_lifespan project
# Batch correction
# Anthony Ramnauth, Apr 21 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(scater)
    library(scran)
    library(harmony)
    library(ggplot2)
    library(ggnewscale)
    library(PCAtools)
    library(spatialLIBD)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "batch_correction")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <-
    readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

# Load BayesSpace k=8 clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "k10_clustering_results"),
    prefix = ""
)

spe$bayesSpace_harmony_10 <- as.factor(spe$bayesSpace_harmony_10)

# Feature selection
dec <- modelGeneVar(spe, block = spe$sample_id)
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
spe <- runPCA(spe, subset_row = top_hvgs, ncomponents = 50)
set.seed(12345)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")

# Plot UMAP before harmony
pdf(file = here::here("plots", "batch_correction", "DG_UMAP_spe.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(spe$sample_id),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "Capture Area") +
    theme_bw() +
    theme_classic()

dev.off()

# run Harmony on PCA dimensions to integrate sample IDs

pca_matrix <- reducedDim(spe, "PCA")
sample_ids <- colData(spe)$sample_id
stopifnot(nrow(pca_matrix) == length(sample_ids))

set.seed(123)
harmony_embeddings <- HarmonyMatrix(pca_matrix,
    meta_data = sample_ids,
    do_pca = FALSE)

colnames(harmony_embeddings) <-
    paste0("HARMONY", seq_len(ncol(harmony_embeddings)))

dim(harmony_embeddings)
head(harmony_embeddings, 2)

# store in SPE object
reducedDims(spe) <- list(
    PCA = reducedDim(spe, "PCA"),
    UMAP = reducedDim(spe, "UMAP"),
    HARMONY = harmony_embeddings
)

# Run UMAP on harmony integrated spe object
set.seed(12345)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

# Plot UMAP after harmony
pdf(file = here::here("plots", "batch_correction", "DG_UMAP_harmony.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(spe$sample_id),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "Capture Area") +
    theme_bw() +
    theme_classic()

dev.off()

# Create directory for harmony integrated spe
dir_rdata <- here::here("processed-data", "harmony_processed_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

saveRDS(spe,
    file = here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))
