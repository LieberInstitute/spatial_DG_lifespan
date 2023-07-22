###################################
# spatial_DG_lifespan project
# Pseudo-bulking Manual Annotations
# Anthony Ramnauth, Oct 11 2022
###################################

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

## Pseudo-bulk for Manual Annotations
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        ManualAnnotation = spe$ManualAnnotation,
        sample_id = spe$sample_id
    )
)

dim(spe_pseudo)

spe_pseudo$ManualAnnotation <- factor(spe_pseudo$ManualAnnotation)

# find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo,
    group = spe_pseudo$ManualAnnotation)

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
# [1] 14.300 12.500  5.360  4.170  3.510  3.100  2.560  2.400  2.330  1.730  1.560  1.510  1.450  1.320
#[15]  1.310  1.270  1.190  1.140  1.110  1.060  1.060  0.966  0.893  0.870  0.841  0.805  0.796  0.710
#[29]  0.677  0.648  0.641  0.633  0.624  0.570  0.561  0.542  0.536  0.520  0.500  0.499  0.478  0.468
#[43]  0.464  0.463  0.453  0.443  0.437  0.424  0.415  0.412

pca_pseudo <- pca$x[, seq_len(50)]
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

percent.var <- jaffelab::getPcaVars(pca)[seq_len(50)]

chosen.elbow <- findElbowPoint(percent.var)
chosen.elbow

# Elbow plot of PCs & plot Reduced Dimensions
pdf(file = here::here("plots", "pseudobulked", "MA_Elbow_plot_spe.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL_CA1" = "#16FF32", "PCL_CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD",
    "SUB" = "#AA0DFE", "THAL" = "navy")

# Plot PCA
pdf(file = here::here("plots", "pseudobulked", "pseudobulked_Manual_Annotated_PCA.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 6, point_size = 1)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 6, point_size = 1)
plotPCA(spe_pseudo, colour_by = "ManualAnnotation", ncomponents = 6, point_size = 1) +
    scale_color_manual(values = man_colors) +
    labs(color = "ManualAnnotation")
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 6, point_size = 1)
dev.off()

#### plot explanatory variables ####

# uses linear regression model
vars <- getVarianceExplained(spe_pseudo,
    variables = c("sex", "age_bin", "ManualAnnotation", "sample_id")
)
head(vars)
#                        sex   age_bin ManualAnnotation sample_id
#ENSG00000241860 0.181785492 0.9022613         16.23925  14.20748
#ENSG00000237491 0.001490725 1.0537814         15.31700  16.94108
#ENSG00000228794 1.717717121 1.2073105         32.97109  18.01712
#ENSG00000225880 2.838839650 1.0993660         11.24023  10.68455
#ENSG00000230368 2.722072681 0.6457584          5.49325  28.71624
#ENSG00000187634 2.467870514 2.7736111         29.48780  15.19095

pdf(file = here::here("plots", "pseudobulked", "plot__manual_annotated_explanatory_vars.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
saveRDS(spe_pseudo, file = here::here("processed-data", "pseudobulk_spe",
    "manual_annotated_pseudobulk_spe.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
