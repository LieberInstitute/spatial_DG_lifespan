#############################################
# spatial_DG_lifespan project
# t-stat correlation of MA & BayesSpace
# Anthony Ramnauth, July 22 2023
#############################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(dplyr)
    library(ComplexHeatmap)
    library(sessioninfo)
})

# Load Maunual Annotation modeling results
MA_modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe",
    "manual_annotated_modeling_results.rds"))

# Load BayesSpace with CP cluster modeling results
BS_modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe",
    "modeling_results_wCP.rds"))

## extract t-statics and rename for BayesSpace enrichment model
registration_t_stats <- BS_modeling_results$enrichment[, grep("^t_stat",
    colnames(BS_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
rownames(registration_t_stats) <- BS_modeling_results$enrichment$ensembl

## spatial domain x gene
dim(registration_t_stats)
# [1] 13355    10

## check out table
registration_t_stats[1:5, 1:5]

cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = MA_modeling_results,
    model_type = "enrichment",
    top_n = 100
)

layer_stat_cor_plot(cor_layer, max = max(cor_layer))

man_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "SO" = "#F6222E", "SR" = "#FE00FA",
    "PCL_CA1" = "#16FF32", "PCL_CA3" = "#3283FE", "CA4" = "#FEAF16", "GCL" = "#B00068",
    "SGZ" = "#1CFFCE", "SL" = "#90AD1C", "WM" = "#2ED9FF", "CP" = "#DEA0FD",
    "SUB" = "#AA0DFE", "THAL" = "navy")

bay_colors <- c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")

pdf(file = here::here("plots", "pseudobulked", "MA_Bayes_tstat_corr.pdf"))

Heatmap(cor_layer,
    name = "cor",
    top_annotation = HeatmapAnnotation(Manual_Anno = colnames(cor_layer),
    col = list(Manual_Anno = man_colors)),
    left_annotation = rowAnnotation(BayesSpace = rownames(cor_layer),
        col = list(BayesSpace = bay_colors)),
    show_column_names = TRUE,
    column_names_rot = 45,
    show_row_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE
    )

dev.off()

