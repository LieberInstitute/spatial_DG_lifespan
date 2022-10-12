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

## Pseudo-bulk for Manual Annotations
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        ManualAnnotation = spe$ManualAnnotation,
        sample_id = spe$sample_id
    )
)

spe_pseudo$ManualAnnotation <- factor(spe_pseudo$ManualAnnotation)

# find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$ManualAnnotation)

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
#[1] 13.800 11.700  6.070  5.040  4.310  3.870  3.310  3.030  2.680  2.300  2.150  2.010  1.900  1.800  1.710  1.610
#[17]  1.550  1.430  1.380  1.290  1.210  1.120  1.050  0.967  0.901  0.891  0.851  0.743  0.706  0.685  0.667  0.664
#[33]  0.617  0.613  0.610  0.592  0.581  0.565  0.555  0.541  0.522  0.487  0.479  0.454  0.442  0.432  0.410  0.402
#[49]  0.401  0.393

pca_pseudo <- pca$x[, seq_len(50)]
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

# Plot PCA
pdf(file = here::here("plots", "pseudobulked", "pseudobulked_Manual_Annotated_PCA.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo, colour_by = "age_bin", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo, colour_by = "ManualAnnotation", ncomponents = 12, point_size = 1)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1)
dev.off()

#### plot explanatory variables ####

# uses linear regression model
vars <- getVarianceExplained(spe_pseudo,
    variables = c("sex", "age_bin", "ManualAnnotation", "sample_id")
)
head(vars)
#                      sex    age_bin ManualAnnotation sample_id
#ENSG00000237491 0.3670081  2.7720179         23.41876  9.382267
#ENSG00000228794 2.2702344  7.2479361         47.10084 16.217047
#ENSG00000225880 3.7682041  0.6288228         26.50485  7.568232
#ENSG00000187634 7.9641812 10.5727159         35.02832 15.797241
#ENSG00000188976 4.9468022 21.0766644         40.48412 25.275481
#ENSG00000187961 2.6804164  5.8721646         30.36934 15.001371

pdf(file = here::here("plots", "pseudobulked", "plot__manual_annotated_explanatory_vars.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
saveRDS(spe_pseudo, file = here::here("processed-data", "pseudobulk_spe", "manual_annotated_pseudobulk_spe.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
