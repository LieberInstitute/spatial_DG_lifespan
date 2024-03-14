#############################################################################
# spatial_DG_lifespan project
# Pseudo-bulking BayesSpace clusters removing low spot count & CP
# Anthony Ramnauth, April 10 2023
#############################################################################

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

dim(spe)

# Remove Visium spots contaminated with Thalamus using
# gene marker TCF7L2 ENSG00000148737
spe <- spe[, which(logcounts(spe)["ENSG00000148737",] <= 1)]

# Remove Choroid Plexus cluster
spe = spe[, which(spe$bayesSpace_harmony_10 != "3")]

# Refactor bayesSpace_harmony_10
spe$bayesSpace_harmony_10 <- as.numeric(spe$bayesSpace_harmony_10)
spe$bayesSpace_harmony_10 <- as.factor(spe$bayesSpace_harmony_10)

dim(spe)

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

dim(spe_pseudo)

# Remove speudo-bulked entries with very low spot counts
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]

dim(spe_pseudo)

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
# [1] 18.200  7.680  5.950  5.200  3.520  3.340  2.740  2.430  2.330  2.240
#[11]  1.840  1.690  1.620  1.570  1.270  1.240  1.220  1.110  1.080  1.000
#[21]  0.983  0.922  0.892  0.881  0.832  0.795  0.764  0.738  0.709  0.701
#[31]  0.670  0.649  0.634  0.619  0.606  0.582  0.562  0.549  0.541  0.522
#[41]  0.513  0.501  0.479  0.471  0.468  0.458  0.437  0.431  0.422  0.402

pca_pseudo <- pca$x[, seq_len(50)]
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

percent.var <- jaffelab::getPcaVars(pca)[seq_len(50)]

chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow

# Elbow plot of PCs & plot Reduced Dimensions
pdf(file = here::here("plots", "pseudobulked", "Elbow_plot_spe.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

# Plot PCA
pdf(file = here::here("plots", "pseudobulked", "pseudobulked_PCA.pdf"))
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

pdf(file = here::here("plots", "pseudobulked", "pseudobulked_PCA2vs1.pdf"),
    width = 12, height = 8)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 10)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 2, point_size = 10) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20))
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 10) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20))
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 2, point_size = 10,
        shape = "age_bin") +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20)) +
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
#                       sex    age_bin BayesSpace sample_id
#ENSG00000237491 0.02079956  2.1107390  14.686445  19.48974
#ENSG00000228794 2.20406900  0.6282807  31.181354  15.53606
#ENSG00000188976 1.93457937  0.9461802  46.606339  20.85111
#ENSG00000187961 0.51088432  3.0575454   8.298991  13.94726
#ENSG00000188290 0.53714589 17.1297162  41.520782  42.90282
#ENSG00000187608 0.20844463 11.3509890  36.773635  50.02008

pdf(file = here::here("plots", "pseudobulked", "plot_explanatory_vars.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
saveRDS(spe_pseudo, file = here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))
