####################################
# spatial_DG_lifespan project
# Pseudo-bulking BayesSpace clusters
# Anthony Ramnauth, May 10 2022
####################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

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
})

# Create directory for BayesSpace pseudo-bulked spe object
dir_rdata <- here::here("processed-data", "pseudobulk_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

# Create directory for pseudobulked plots
dir_plots <- here::here("plots", "pseudobulked")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Remove Visium spots contaminated with Thalamus using
# gene marker TCF7L2 ENSG00000148737
spe <- spe[, which(logcounts(spe)["ENSG00000148737",] <= 1)]

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        between(spe.age, 0, 3) ~ "Infant",
        between(spe.age, 13, 19) ~ "Teen",
        between(spe.age, 20, 50) ~ "Adult",
        between(spe.age, 60, 100) ~ "Elderly"
    ))

stopifnot(age_df$spe.key == spe$key)

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

## Pseudo-bulk for BayesSpace k = 10 results
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        BayesSpace = spe$bayesSpace_harmony_10,
        sample_id = spe$sample_id
    )
)

# find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$BayesSpace)

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
# [1] 16.800 13.000  5.060  4.480  3.200  2.840  2.760  2.660  2.530  2.340
#[11]  2.030  1.980  1.720  1.550  1.480  1.390  1.320  1.250  1.180  1.140
#[21]  1.090  0.972  0.928  0.913  0.866  0.830  0.798  0.763  0.681  0.653
#[31]  0.646  0.572  0.565  0.525  0.516  0.510  0.477  0.460  0.460  0.440
#[41]  0.421  0.402  0.390  0.380  0.372  0.356  0.346  0.336  0.334  0.323

pca_pseudo <- pca$x[, seq_len(50)]
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

percent.var <- jaffelab::getPcaVars(pca)[seq_len(50)]

chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow

# Elbow plot of PCs & plot Reduced Dimensions
pdf(file = here::here("plots", "pseudobulked", "Elbow_plot_spe_wCP.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

# Plot PCA
pdf(file = here::here("plots", "pseudobulked", "pseudobulked_PCA_wLowSpots.pdf"))
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "race", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "rin", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "pmi", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 6) +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace")
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 6)
plotPCA(spe_pseudo, colour_by = "ncells", ncomponents = 6)
dev.off()

pdf(file = here::here("plots", "pseudobulked", "pseudobulked_PCA2vs1_wLowSpots.pdf"))
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 10)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 2, point_size = 10)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 10)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 2, point_size = 10) +
    scale_color_manual(values = bay_colors) +
    labs(color = "BayesSpace")
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 2, point_size = 10)
dev.off()

#### plot explanatory variables ####

# uses linear regression model
vars <- getVarianceExplained(spe_pseudo,
    variables = c("sex", "age_bin", "BayesSpace", "sample_id")
)
head(vars)
#                        sex   age_bin BayesSpace sample_id
#ENSG00000237491 0.513646905 1.2581478   12.21868 15.254694
#ENSG00000228794 0.617999458 0.8668777   26.43471 12.073696
#ENSG00000188976 2.072129847 1.6700312   25.88816 19.182095
#ENSG00000187961 1.128175139 1.6840020   17.75967  9.300544
#ENSG00000188290 0.081362291 8.4818055   21.00527 28.188464
#ENSG00000187608 0.004644758 3.8556762   25.96791 39.477292

pdf(file = here::here("plots", "pseudobulked", "plot_explanatory_vars_wLowSpots.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
saveRDS(spe_pseudo, file = here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe_wLowSpots.rds"))
