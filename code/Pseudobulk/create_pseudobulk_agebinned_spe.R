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
spe <- readRDS(here::here("processed-data", "BayesSpace", "bayesspace_first_two_slides_spe.rds"))

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
sce <- as(spe, "SingleCellExperiment")
spe_pseudo <- aggregateAcrossCells(
    sce,
    DataFrame(
        BayesSpace = spe$BayesSpace8,
        age_bin = spe$age_bin
    )
)


# find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_age_bin <- filterByExpr(spe_pseudo, group = spe_pseudo$age_bin)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$BayesSpace)

summary(rowData(spe_pseudo)$high_expr)
summary(rowData(spe_pseudo)$high_expr_group_age_bin)
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
#[1] 3.12e+01 1.13e+01 8.59e+00 8.29e+00 6.56e+00 5.47e+00 3.73e+00 3.37e+00 2.83e+00 2.14e+00 1.88e+00 1.76e+00
#[13] 1.37e+00 1.13e+00 1.00e+00 9.81e-01 9.59e-01 8.61e-01 7.43e-01 6.99e-01 6.52e-01 5.92e-01 5.71e-01 5.49e-01
#[25] 5.08e-01 4.67e-01 4.54e-01 3.99e-01 3.18e-01 3.04e-01 2.69e-01 5.39e-29       NA       NA       NA       NA
#[37]       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA       NA
#[49]       NA       NA

pca_pseudo <- pca$x[, seq_len(32)]
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

# Plot PCA
pdf(file = here::here("plots", "pseudobulked", "First_two_slides_pseudobulked_PCA.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 15, point_size = 1)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 15, point_size = 1)
plotPCA(spe_pseudo, colour_by = "ncells", ncomponents = 15, point_size = 1)
dev.off()

#### plot explanatory variables ####

# uses linear regression model
vars <- getVarianceExplained(spe_pseudo,
    variables = c("age_bin", "BayesSpace", "ncells")
)
head(vars)
#                  age_bin BayesSpace      ncells
#ENSG00000237491 12.783255   36.11698 14.66670805
#ENSG00000228794  4.925178   63.20400  3.65738242
#ENSG00000187634 24.844835   52.22054  0.01972499
#ENSG00000188976 19.303919   68.68344  0.33715353
#ENSG00000187961 23.173932   19.51365  0.43828083
#ENSG00000188290 12.391068   75.71368  0.04083123

pdf(file = here::here("plots", "pseudobulked", "First_second_slide_plot_explanatory_vars.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
saveRDS(spe_pseudo, file = here::here("processed-data", "pseudobulk_spe", "First_two_pseudobulk_spe.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
