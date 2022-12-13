####################################
# spatial_DG_lifespan project
# Pseudo-bulking BayesSpace clusters
# Anthony Ramnauth, May 10 2022
####################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(edgeR)
    library(scuttle)
    library(scater)
    library(scran)
    library(dplyr)
    library(PCAtools)
    library(sessioninfo)
})

# Create directory for BayesSpace pseudo-bulked spe object
dir_rdata <- here::here("processed-data", "pseudobulk_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

# Create directory for pseudobulked plots
dir_plots <- here::here("plots", "pseudobulked")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

# Remove CP cluster
spe = spe[, which(spe$bayesSpace_harmony_8 != "5")]

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        grepl("Br1412", age_df$spe.sample_id) ~ "Teen",
        grepl("Br2706", age_df$spe.sample_id) ~ "Teen",
        grepl("Br3942", age_df$spe.sample_id) ~ "Adult",
        grepl("Br5242", age_df$spe.sample_id) ~ "Elderly",
        grepl("Br6023", age_df$spe.sample_id) ~ "Elderly",
        grepl("Br8195", age_df$spe.sample_id) ~ "Infant",
        grepl("Br8667", age_df$spe.sample_id) ~ "Adult",
        grepl("Br8686", age_df$spe.sample_id) ~ "Infant"
    ))

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

## Pseudo-bulk for BayesSpace k = 8 results
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        BayesSpace = spe$bayesSpace_harmony_8,
        sample_id = spe$sample_id
    )
)

spe_pseudo$BayesSpace <- factor(spe_pseudo$BayesSpace)

# find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$bayesSpace_harmony_8)

summary(rowData(spe_pseudo)$high_expr)
summary(rowData(spe_pseudo)$high_expr_group_sample_id)
summary(rowData(spe_pseudo)$high_expr_group_cluster)

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)

# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x

dim(spe_pseudo)

# run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))
jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 22.900 11.000  7.620  6.160  5.470  4.410  4.270  3.780  3.030  2.620  2.140  1.940  1.720  1.310  1.200  1.130
#[17]  1.040  1.020  0.982  0.887  0.851  0.806  0.779  0.719  0.692  0.673  0.625  0.621  0.612  0.579  0.555  0.534
#[33]  0.513  0.499  0.474  0.450  0.435  0.414  0.398  0.359  0.341  0.327  0.315  0.304  0.287  0.276  0.274  0.268
#[49]  0.248  0.236

pca_pseudo <- pca$x[, seq_len(50)]
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

# Plot PCA
pdf(file = here::here("plots", "pseudobulked", "pseudobulked_PCA.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1)
dev.off()

#### plot explanatory variables ####

# uses linear regression model
vars <- getVarianceExplained(spe_pseudo,
    variables = c("sex", "age_bin", "BayesSpace", "sample_id")
)
head(vars)
#                      sex   age_bin BayesSpace sample_id
#ENSG00000237491 0.7277197  7.969127   22.40622  17.24963
#ENSG00000228794 1.8585308  6.314668   47.31703  10.35570
#ENSG00000188976 4.0273672 15.993802   48.20171  20.17645
#ENSG00000187961 5.1299084  2.817930   12.39323  20.09976
#ENSG00000188290 1.1317946  2.884313   46.67667  38.11036
#ENSG00000187608 1.9757806 19.671871   29.13919  59.43019

pdf(file = here::here("plots", "pseudobulked", "plot_explanatory_vars.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
saveRDS(spe_pseudo, file = here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
