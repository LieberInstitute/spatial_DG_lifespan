################################
# spatial_DG_lifespan project
# Batch correction of sce object
# Anthony Ramnauth, Nov 05 2022
################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(scater)
    library(scran)
    library(harmony)
    library(ggplot2)
    library(ggnewscale)
    library(PCAtools)
    library(spatialLIBD)
    library(sessioninfo)
})

# Load SCE
sce <- readRDS(here::here("processed-data", "QC_processed_sce", "QCed_sce.rds"))

# Feature selection
dec <- modelGeneVar(sce, block = sce$Dataset)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
head(top.hvgs.fdr5)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
head(top.hvgs.fdr1)

# Dimensionality reduction
set.seed(12345)
sce <- runPCA(sce, subset_row = top_hvgs, ncomponents = 50)
set.seed(12345)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- c("UMAP1", "UMAP2")

# Plot UMAP before harmony
pdf(file = here::here("sce_plots", "DG_UMAP_sce.pdf"))
ggplot(
    data.frame(reducedDim(sce, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(sce$Dataset), alpha = 0.01)
) +
    geom_point() +
    labs(color = "Dataset") +
    theme_bw()

dev.off()

# run Harmony on PCA dimensions to integrate sample IDs

pca_matrix <- reducedDim(sce, "PCA")
Dataset <- colData(sce)$Dataset
stopifnot(nrow(pca_matrix) == length(Dataset))

set.seed(123)
harmony_embeddings <- HarmonyMatrix(
    pca_matrix,
    meta_data = Dataset,
    do_pca = FALSE
)

colnames(harmony_embeddings) <- paste0("HARMONY", seq_len(ncol(harmony_embeddings)))

dim(harmony_embeddings)
head(harmony_embeddings, 2)

# store in SCE object
reducedDims(sce) <- list(
    PCA = reducedDim(sce, "PCA"),
    UMAP = reducedDim(sce, "UMAP"),
    HARMONY = harmony_embeddings
)

# Run UMAP on harmony integrated sce object
set.seed(123)
sce <- runUMAP(sce, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

# Plot UMAP after harmony
pdf(file = here::here("sce_plots", "DG_UMAP_harmony_sce.pdf"))
ggplot(
    data.frame(reducedDim(sce, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(sce$Dataset), alpha = 0.01)
) +
    geom_point() +
    labs(color = "Dataset") +
    theme_bw()

dev.off()

# Save new sce object
saveRDS(sce, file = here::here("sce_objects", "sce_harmony.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
