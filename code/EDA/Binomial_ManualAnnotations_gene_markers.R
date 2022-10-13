#########################################
# spatial_DG_lifespan project
# Top Marker Genes for Manual Annotations
# Anthony Ramnauth, Oct 10 2022
#########################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(spatialLIBD)
    library(SpatialExperiment)
    library(BayesSpace)
    library(scran)
    library(scater)
    library(dplyr)
    library(pheatmap)
    library(ComplexHeatmap)
    library(circlize)
})

# Create directory for Top Marker Genes plots
dir_plots <- here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

#################################################
# Test for marker genes across all tissue samples
#################################################

markers <- findMarkers(spe, groups = spe$ManualAnnotation, test = "binom", direction = "up")

# Returns a list with one DataFrame per Manual Annotation
markers

###########################################################
# Make .csv lists of top markers for each Manual Annotation
###########################################################

# Make a data frame summary
binom_GCL <- markers[["GCL"]]
GCL_binom_summary <- data.frame(
    gene_name = rownames(binom_GCL),
    rank = binom_GCL$Top,
    p_val = binom_GCL$p.value,
    FDR = binom_GCL$FDR
)

GCL_binom_summary <- GCL_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
dir_outputs <- here("processed-data", "spatialLIBD_manual_annotations")
fn_out_1 <- file.path(dir_outputs, "GCL_binomial_test_results")

# Export summary as .csv file
write.csv(GCL_binom_summary,fn_out_1, row.names = FALSE)

# Make a data frame summary
binom_ML <- markers[["ML"]]
ML_binom_summary <- data.frame(
    gene_name = rownames(binom_ML),
    rank = binom_ML$Top,
    p_val = binom_ML$p.value,
    FDR = binom_ML$FDR
)

ML_binom_summary <- ML_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_2 <- file.path(dir_outputs, "ML_binomial_test_results")

# Export summary as .csv file
write.csv(ML_binom_summary,fn_out_2, row.names = FALSE)

# Make a data frame summary
binom_SGZ <- markers[["SGZ"]]
SGZ_binom_summary <- data.frame(
    gene_name = rownames(binom_SGZ),
    rank = binom_SGZ$Top,
    p_val = binom_SGZ$p.value,
    FDR = binom_SGZ$FDR
)

SGZ_binom_summary <- SGZ_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_3 <- file.path(dir_outputs, "SGZ_binomial_test_results")

# Export summary as .csv file
write.csv(SGZ_binom_summary,fn_out_3, row.names = FALSE)

# Make a data frame summary
binom_CA4 <- markers[["CA4"]]
CA4_binom_summary <- data.frame(
    gene_name = rownames(binom_CA4),
    rank = binom_CA4$Top,
    p_val = binom_CA4$p.value,
    FDR = binom_CA4$FDR
)

CA4_binom_summary <- CA4_binom_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_CA4 <- file.path(dir_outputs, "CA4_binomial_test_results")

# Export summary as .csv file
write.csv(CA4_binom_summary,fn_out_CA4, row.names = FALSE)

###############################################################
# Plot log-fold changes for one cluster over all other clusters
###############################################################

f1 = colorRamp2(c(-5, 0, 5), c("blue", "#EEEEEE", "red"))

# Plot heatmap of GCL top markers
interesting_GCL <- markers[["GCL"]]
best_set_GCL <- interesting_GCL[interesting_GCL$Top <= 5, ]
logFCs_GCL <- getMarkerEffects(best_set_GCL)

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "binomial_for_GCL.pdf"))
Heatmap(logFCs_GCL,
    name = "logFC",
    col = f1,
    column_title = "logFC for GCL from binomial test")
dev.off()

# Plot log-transformed normalized expression of top genes for GCL
top_genes_GCL <- head(rownames(interesting_GCL))

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "top_genes_for_binomial_GCL.pdf"))
plotExpression(spe,
    x = "ManualAnnotation",
    features = top_genes_GCL,
    colour_by = "ManualAnnotation",
    ) +
    theme(axis.text.x=element_blank())
dev.off()

# Plot heatmap of ML top markers
interesting_ML <- markers[["ML"]]
best_set_ML <- interesting_ML[interesting_ML$Top <= 5, ]
logFCs_ML <- getMarkerEffects(best_set_ML)

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "binomial_for_ML.pdf"))
pheatmap(logFCs_ML,
    name = "logFC",
    col = f1,
    main = "logFC for ML from binomial test")
dev.off()

# Plot log-transformed normalized expression of top genes for ML
top_genes_ML <- head(rownames(interesting_ML))

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "top_genes_for_binomial_ML.pdf"))
plotExpression(spe,
    x = "ManualAnnotation",
    features = top_genes_ML,
    colour_by = "ManualAnnotation",
    ) +
    theme(axis.text.x=element_blank())
dev.off()

# Plot heatmap of SGZ top markers
interesting_SGZ <- markers[["SGZ"]]
best_set_SGZ <- interesting_SGZ[interesting_SGZ$Top <= 5, ]
logFCs_SGZ <- getMarkerEffects(best_set_SGZ)

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "binomial_for_SGZ.pdf"))
pheatmap(logFCs_SGZ,
    name = "logFC",
    col = f1,
    main = "logFC for SGZ from binomial test")
dev.off()

# Plot log-transformed normalized expression of top genes for SGZ
top_genes_SGZ <- head(rownames(interesting_SGZ))

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "top_genes_for_binomial_SGZ.pdf"))
plotExpression(spe,
    x = "ManualAnnotation",
    features = top_genes_SGZ,
    colour_by = "ManualAnnotation",
    ) +
    theme(axis.text.x=element_blank())
dev.off()

# Plot heatmap of CA4 top markers
interesting_CA4 <- markers[["CA4"]]
best_set_CA4 <- interesting_CA4[interesting_CA4$Top <= 5, ]
logFCs_CA4 <- getMarkerEffects(best_set_CA4)

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "binomial_for_CA4.pdf"))
pheatmap(logFCs_CA4,
    name = "logFC",
    col = f1,
    main = "logFC for CA4 from binomial test")
dev.off()

# Plot log-transformed normalized expression of top genes for CA4
top_genes_CA4 <- head(rownames(interesting_CA4))

pdf(file = here::here("plots", "manual_annotations", "Binomial_ManualAnnotations_genes",
    "top_genes_for_binomial_CA4.pdf"))
plotExpression(spe,
    x = "ManualAnnotation",
    features = top_genes_CA4,
    colour_by = "ManualAnnotation",
    ) +
    theme(axis.text.x=element_blank())
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
