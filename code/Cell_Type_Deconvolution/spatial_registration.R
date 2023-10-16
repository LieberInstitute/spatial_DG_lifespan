######################################################################
# spatial_DG_lifespan project
# Spatial Registration of formatted Franjic et al., 2021 snRNASeq data
# Anthony Ramnauth, Sept 28 2023
######################################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(spatialLIBD)
    library(SingleCellExperiment)
    library(dplyr)
})

# Load SCE
sce <- readRDS(file = here::here("processed-data", "sce", "sce_sestan_DG_final.rds"))

# Load modeling results from lifespan DG Visium
modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

gene_ensembl <- rownames(sce)
rowData(sce)$gene_ensembl <- gene_ensembl

## Perform the spatial registration
sce_modeling_results <- registration_wrapper(
    sce = sce,
    var_registration = "Cell_Type",
    var_sample_id = "sample_ID",
    gene_ensembl = "gene_ensembl",
    gene_name = "gene_name"
)

## extract t-statics and rename
registration_t_stats <- sce_modeling_results$enrichment[, grep("^t_stat", colnames(sce_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

## cell types x gene
dim(registration_t_stats)
# 19290    22

## check out table
registration_t_stats[1:5, 1:5]
#                  Astro_1     Astro_2      CA2_N        CA3_N        COP
# ENSG00000238009 1.399124 -0.09874905  0.4839145  1.736320260 -0.3165016
# ENSG00000237491 0.315798  0.46802508  0.6352494  1.693940483 -0.2485864
# ENSG00000225880 1.370195  1.72110165 -1.1437533 -0.616894966 -2.3661126
# ENSG00000272438 4.445627  1.55036765  1.2462646 -0.005029478 -1.3940950
# ENSG00000187634 8.507269  0.99030243 -0.7099496 -1.341592790 -1.5365877

cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = modeling_results,
    model_type = "enrichment",
    top_n = 100
)

cor_layer
#                   CA1       CA3_4         GCL          ML          SGZ           SL         SLM          SR         WM
# Astro_1   -0.38767062 -0.39845751 -0.21806843  0.36951781 -0.205011160 -0.094859031  0.55966949  0.36206932  0.2091726
# Astro_2   -0.35747936 -0.36362085 -0.23121553  0.40793335 -0.219576795 -0.048987494  0.51663905  0.37144960  0.1638272
# Pericyte  -0.12979861 -0.24042492 -0.31189314  0.07639730 -0.287897586  0.236476448  0.32812865  0.38565041  0.1948331
# Endoth    -0.08173001 -0.20400719 -0.21222217  0.08824174 -0.229334188  0.155339569  0.24699464  0.31353867  0.1187624
# SMC       -0.10295299 -0.22171027 -0.27172615  0.05439505 -0.185192278  0.144382282  0.28750430  0.33603936  0.1732441
# Oligo     -0.19827787 -0.18482811 -0.08219674 -0.21762882 -0.376168104 -0.005772486  0.16723343  0.04095009  0.6723950
# VLMC      -0.20514492 -0.20317116 -0.23087092  0.06382088 -0.112488651  0.001731823  0.29250002  0.31155412  0.2242552
# Macro     -0.17019814 -0.17664560 -0.10173192  0.04233869  0.004082427 -0.049083018  0.18173182  0.19152388  0.1581640
# Microglia -0.14931665 -0.19457393 -0.12994998  0.08455658 -0.125496609  0.064200130  0.20195840  0.22143698  0.1488086
# T_Cell    -0.04111028 -0.05383486 -0.13665020 -0.04701647 -0.155489494  0.154705674  0.09728729  0.13758983  0.1209615
# COP       -0.12229694 -0.11765458 -0.09861969 -0.06104832 -0.046608996 -0.021745194  0.13817259  0.08413711  0.2557846
# OPC       -0.11635187 -0.13441360 -0.01567185  0.09709994 -0.022628511 -0.062912833  0.11545429  0.05410720  0.1007566
# GC         0.12540403  0.16918355  0.77968516  0.36297644 -0.139850794 -0.361980959 -0.46690973 -0.51392407 -0.4237145
# CA3_N      0.48463448  0.70558028  0.31992477 -0.19487618  0.013704048 -0.003373391 -0.64454345 -0.51661636 -0.4866116
# CA2_N      0.38028563  0.59130794  0.28912502 -0.16545434  0.047334952 -0.053616619 -0.54205975 -0.44911864 -0.3995773
# Mossy      0.36074142  0.60815237  0.28836387 -0.18809767  0.180325352 -0.120583235 -0.56363209 -0.47266648 -0.4021644
# InN_SST    0.26014470  0.23699759  0.07276447 -0.30650754  0.637232424  0.053856154 -0.36091812 -0.31130603 -0.2751321
# InN_VIP    0.19135281  0.23032509  0.20753627 -0.12652957  0.427123051 -0.077431701 -0.32634990 -0.28817001 -0.3219581
# InN_LAMP5  0.20436058  0.17523920  0.07836634 -0.19066276  0.477103778 -0.015294109 -0.23834594 -0.19728974 -0.2745286
# InN_LHX6   0.27231217  0.21832298  0.13173421 -0.23310883  0.371306685  0.036746256 -0.32665966 -0.28286785 -0.2462427
# InN_NR2F2  0.24699449  0.26104230  0.12248242 -0.23694771  0.414529309  0.031864275 -0.31586678 -0.29132939 -0.2888577
# InN_PV     0.25655840  0.23216851  0.11833563 -0.20695281  0.431044294  0.028822947 -0.33395540 -0.26436578 -0.2846175


pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Spatial_registration.pdf"))

layer_stat_cor_plot(cor_layer, max = max(cor_layer))

dev.off()

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno

# Scatterplots of t-stats for snRNAseq cell type vs Visuim spatial domain

# Exclude ensemblIDs of snRNAseq not in Visium dataset and vice versa

v_registration_t_stats <- registration_t_stats[rownames(registration_t_stats) %in% rownames(modeling_results$enrichment),]
v_modeling_results <- modeling_results$enrichment[rownames(modeling_results$enrichment) %in% rownames(registration_t_stats),]

dim(v_registration_t_stats)
dim(v_modeling_results)

#Reorder rows so ensemblIDs match
v_registration_t_stats <- v_registration_t_stats[! v_registration_t_stats %in%
        setdiff(rownames(v_registration_t_stats), rownames(v_modeling_results)),]

v_registration_t_stats <- v_registration_t_stats[order(match(rownames(v_registration_t_stats), rownames(v_modeling_results))), ]

stopifnot(rownames(v_registration_t_stats) == rownames(v_modeling_results))

# Plot

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "tstat_scatterplot.pdf"))

plot(v_modeling_results$t_stat_GCL, v_registration_t_stats$GC)
abline(lm(v_registration_t_stats$GC ~ v_modeling_results$t_stat_GCL), col = "blue")

plot(v_modeling_results$t_stat_WM, v_registration_t_stats$Oligo)
abline(lm(v_registration_t_stats$Oligo ~ v_modeling_results$t_stat_WM), col = "blue")

plot(v_modeling_results$t_stat_CA3_4, v_registration_t_stats$CA3_N)
abline(lm(v_registration_t_stats$CA3_N ~ v_modeling_results$t_stat_CA3_4), col = "blue")

plot(v_modeling_results$t_stat_CA3_4, v_registration_t_stats$Mossy)
abline(lm(v_registration_t_stats$Mossy ~ v_modeling_results$t_stat_CA3_4), col = "blue")

plot(v_modeling_results$t_stat_SGZ, v_registration_t_stats$InN_SST)
abline(lm(v_registration_t_stats$InN_SST ~ v_modeling_results$t_stat_SGZ), col = "blue")

plot(v_modeling_results$t_stat_SGZ, v_registration_t_stats$InN_VIP)
abline(lm(v_registration_t_stats$InN_VIP ~ v_modeling_results$t_stat_SGZ), col = "blue")

plot(v_modeling_results$t_stat_SGZ, v_registration_t_stats$InN_LAMP5)
abline(lm(v_registration_t_stats$InN_LAMP5 ~ v_modeling_results$t_stat_SGZ), col = "blue")

plot(v_modeling_results$t_stat_SLM, v_registration_t_stats$Astro_1)
abline(lm(v_registration_t_stats$Astro_1 ~ v_modeling_results$t_stat_SLM), col = "blue")

plot(v_modeling_results$t_stat_SLM, v_registration_t_stats$Astro_2)
abline(lm(v_registration_t_stats$Astro_2 ~ v_modeling_results$t_stat_SLM), col = "blue")

dev.off()
