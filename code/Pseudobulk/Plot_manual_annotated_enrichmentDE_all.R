###############################
# spatial_DG_lifespan project
# Plot DE of Manual Annotations
# Anthony Ramnauth, Oct 12 2022
###############################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(pheatmap)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "manual_annotated_pseudobulk_spe.rds"))

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe",
    "manual_annotated_modeling_results.rds"))

################################################################
# Make a .csv list of enriched genes for each Manual Annotations
################################################################

# Make a data frame summary
CA4_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_CA4,
    p_val = modeling_results$enrichment$p_value_CA4,
    FDR = modeling_results$enrichment$fdr_CA4
)

CA4_enriched_summary <- CA4_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

CA4_enriched_summary <- arrange(CA4_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "spatialLIBD_manual_annotations")
fn_out_CA4 <- file.path(dir_outputs, "CA4_enriched_results")

# Export summary as .csv file
write.csv(CA4_enriched_summary,fn_out_CA4, row.names = FALSE)

# Make a data frame summary
CP_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_CP,
    p_val = modeling_results$enrichment$p_value_CP,
    FDR = modeling_results$enrichment$fdr_CP
)

CP_enriched_summary <- CP_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

CP_enriched_summary <- arrange(CP_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_CP <- file.path(dir_outputs, "CP_enriched_results")

# Export summary as .csv file
write.csv(CP_enriched_summary,fn_out_CP, row.names = FALSE)

# Make a data frame summary
GCL_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_GCL,
    p_val = modeling_results$enrichment$p_value_GCL,
    FDR = modeling_results$enrichment$fdr_GCL
)

GCL_enriched_summary <- GCL_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

GCL_enriched_summary <- arrange(GCL_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_GCL<- file.path(dir_outputs, "GCL_enriched_results")

# Export summary as .csv file
write.csv(GCL_enriched_summary,fn_out_GCL, row.names = FALSE)

# Make a data frame summary
ML_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_ML,
    p_val = modeling_results$enrichment$p_value_ML,
    FDR = modeling_results$enrichment$fdr_ML
)

ML_enriched_summary <- ML_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

ML_enriched_summary <- arrange(ML_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_ML<- file.path(dir_outputs, "ML_enriched_results")

# Export summary as .csv file
write.csv(ML_enriched_summary,fn_out_ML, row.names = FALSE)

# Make a data frame summary
CA1_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_PCL_CA1,
    p_val = modeling_results$enrichment$p_value_PCL_CA1,
    FDR = modeling_results$enrichment$fdr_PCL_CA1
)

CA1_enriched_summary <- CA1_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

CA1_enriched_summary <- arrange(CA1_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_CA1<- file.path(dir_outputs, "CA1_enriched_results")

# Export summary as .csv file
write.csv(CA1_enriched_summary,fn_out_CA1, row.names = FALSE)

# Make a data frame summary
CA3_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_PCL_CA3,
    p_val = modeling_results$enrichment$p_value_PCL_CA3,
    FDR = modeling_results$enrichment$fdr_PCL_CA3
)

CA3_enriched_summary <- CA3_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

CA3_enriched_summary <- arrange(CA3_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_CA3<- file.path(dir_outputs, "CA3_enriched_results")

# Export summary as .csv file
write.csv(CA3_enriched_summary,fn_out_CA3, row.names = FALSE)

# Make a data frame summary
SGZ_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_SGZ,
    p_val = modeling_results$enrichment$p_value_SGZ,
    FDR = modeling_results$enrichment$fdr_SGZ
)

SGZ_enriched_summary <- SGZ_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SGZ_enriched_summary <- arrange(SGZ_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SGZ<- file.path(dir_outputs, "SGZ_enriched_results")

# Export summary as .csv file
write.csv(SGZ_enriched_summary,fn_out_SGZ, row.names = FALSE)

# Make a data frame summary
SL_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_SL,
    p_val = modeling_results$enrichment$p_value_SL,
    FDR = modeling_results$enrichment$fdr_SL
)

SL_enriched_summary <- SL_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SL_enriched_summary <- arrange(SL_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SL<- file.path(dir_outputs, "SL_enriched_results")

# Export summary as .csv file
write.csv(SL_enriched_summary,fn_out_SL, row.names = FALSE)

# Make a data frame summary
SLM_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_SLM,
    p_val = modeling_results$enrichment$p_value_SLM,
    FDR = modeling_results$enrichment$fdr_SLM
)

SLM_enriched_summary <- SLM_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SLM_enriched_summary <- arrange(SLM_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SLM<- file.path(dir_outputs, "SLM_enriched_results")

# Export summary as .csv file
write.csv(SLM_enriched_summary,fn_out_SLM, row.names = FALSE)

# Make a data frame summary
SO_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_SO,
    p_val = modeling_results$enrichment$p_value_SO,
    FDR = modeling_results$enrichment$fdr_SO
)

SO_enriched_summary <- SO_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SO_enriched_summary <- arrange(SO_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SO<- file.path(dir_outputs, "SO_enriched_results")

# Export summary as .csv file
write.csv(SO_enriched_summary,fn_out_SO, row.names = FALSE)

# Make a data frame summary
SR_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_SR,
    p_val = modeling_results$enrichment$p_value_SR,
    FDR = modeling_results$enrichment$fdr_SR
)

SR_enriched_summary <- SR_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SR_enriched_summary <- arrange(SR_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SR<- file.path(dir_outputs, "SR_enriched_results")

# Export summary as .csv file
write.csv(SR_enriched_summary,fn_out_SR, row.names = FALSE)

# Make a data frame summary
SUB_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_SUB,
    p_val = modeling_results$enrichment$p_value_SUB,
    FDR = modeling_results$enrichment$fdr_SUB
)

SUB_enriched_summary <- SUB_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

SUB_enriched_summary <- arrange(SUB_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_SUB<- file.path(dir_outputs, "SUB_enriched_results")

# Export summary as .csv file
write.csv(SUB_enriched_summary,fn_out_SUB, row.names = FALSE)

# Make a data frame summary
WM_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_WM,
    p_val = modeling_results$enrichment$p_value_WM,
    FDR = modeling_results$enrichment$fdr_WM
)

WM_enriched_summary <- WM_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

WM_enriched_summary <- arrange(WM_enriched_summary, desc(t_stat))

# directory to save lists
fn_out_WM<- file.path(dir_outputs, "WM_enriched_results")

# Export summary as .csv file
write.csv(WM_enriched_summary,fn_out_WM, row.names = FALSE)

################################################
# Plot top enriched genes per BayesSpace Cluster
################################################

# Get mean expression
mat <- assays(spe_pseudo)$logcounts

# filter
gIndex <- rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter <- mat[gIndex, ] # subset matrix on just those genes.  want to remove lowly expressed genes.

# Extract the p-values
pvals <- modeling_results$enrichment[, 14:26]
rownames(pvals) <- rownames(mat_filter)

# Extract the t-statistics
t_stat <- modeling_results$enrichment[, 1:13]
rownames(t_stat) <- rownames(mat_filter)

# Extract the FDRs
fdrs <- modeling_results$enrichment[, 27:39]
rownames(fdrs) <- rownames(mat_filter)

### pick top 10 genes per cluster:sample
cluster_specific_indices <- mapply(
    function(t, p, f) {
        oo <- order(t, decreasing = TRUE)[1:6]
    },
    as.data.frame(t_stat),
    as.data.frame(pvals),
    as.data.frame(fdrs)
)
cluster_ind <- unique(as.numeric(cluster_specific_indices))

# Add logcounts from indexed from top genes
exprs_heatmap <- assays(spe_pseudo)[[2]][cluster_ind, ]
rownames(exprs_heatmap) <- rowData(spe_pseudo)$gene_name[cluster_ind]
colnames(exprs_heatmap) <- paste("logcount", 1:83, sep = "")

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "pseudobulked",
    "manual_annotations_enrichment_heatmap_all.pdf"), width = 12, height = 8)
Heatmap(exprs_heatmap,
    name = "logcounts",
    top_annotation = HeatmapAnnotation(Manual_Annotations = spe_pseudo$ManualAnnotation, age = spe_pseudo$age_bin,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen"),
        Manual_Annotations = c("CA4" = "orangered", "CP" = "orange", "GCL" = "cyan", "ML" = "springgreen3",
            "PCL_CA1" = "brown", "PCL_CA3" = "pink", "SGZ" = "yellow", "SL" = "slategrey", "SLM" = "chartreuse",
            "SO" = "coral1", "SR" = "blanchedalmond", "SUB" = "darkseagreen2", "WM" = "black"))),
    column_title = "logcounts of top 6 transcripts from Differential Enrichment of Manual Annotations",
    show_column_names = FALSE,
    row_title = NULL,
    column_split = 18,
    row_split =29,
    row_names_gp = gpar(fontsize = 7)
    )
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
