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
# [1] 22.700  8.810  6.620  5.440  4.980  4.200  3.950  3.370  3.300  2.890  2.650  2.520  2.310  2.040  1.870  1.650
#[17]  1.500  1.400  1.220  1.090  0.891  0.797  0.751  0.678  0.666  0.628  0.601  0.557  0.549  0.525  0.500  0.460
#[33]  0.447  0.437  0.409  0.404  0.392  0.378  0.358  0.342  0.333  0.327  0.310  0.291  0.279  0.270  0.261  0.237
#[49]  0.221  0.215

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
#                       sex   age_bin BayesSpace sample_id
#ENSG00000228794 0.02447535  3.669668   41.22758  7.803956
#ENSG00000188976 1.72967711 12.459186   51.33415 13.380178
#ENSG00000188290 0.02718252  3.095280   41.84710 33.058147
#ENSG00000187608 0.73643273 16.397670   26.07864 50.631397
#ENSG00000188157 3.75328432  4.023631   25.99233 23.633593
#ENSG00000131591 1.17400898 22.869309   19.61353 29.150242

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
