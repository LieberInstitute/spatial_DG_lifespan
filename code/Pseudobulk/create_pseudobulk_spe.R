###############################################
# spatial_DG_lifespan project
# Pseudo-bulking BayesSpace clusters by age_bin
# Anthony Ramnauth, May 10 2022
###############################################

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
# [1] 21.700  8.270  7.580  6.030  5.130  4.770  4.010  3.910  3.380  3.150  2.980
# [12]  2.510  2.130  1.900  1.630  1.470  1.420  1.340  1.150  1.050  0.895  0.793
# [23]  0.713  0.674  0.600  0.587  0.555  0.523  0.496  0.484  0.462  0.419  0.409
# [34]  0.401  0.373  0.372  0.361  0.345  0.338  0.325  0.310  0.297  0.280  0.276
# [45]  0.252  0.244  0.237  0.221  0.204  0.203

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
#                         sex   age_bin BayesSpace sample_id
# ENSG00000228794 2.712823e-01  1.296835  39.083536  6.420347
# ENSG00000188976 1.570477e+00 15.257619  52.724485 17.135663
# ENSG00000188290 4.670369e-04  3.031042  40.019717 30.285038
# ENSG00000187608 6.755845e-01 18.057959  27.046809 52.868439
# ENSG00000188157 2.467874e+00  5.386343  27.617885 21.939271
# ENSG00000131591 1.635477e-03 10.961085  27.675058 20.291421

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
