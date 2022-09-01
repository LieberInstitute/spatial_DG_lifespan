##################################
# spatial_DG_lifespan project
# QC individual Visium Slide
# Anthony Ramnauth, August 22 2022
##################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(scater)
    library(ggplot2)
    library(ggnewscale)
    library(ggspavis)
    library(scuttle)
    library(scran)
    library(dplyr)
    library(sessioninfo)
})

# load saved SPE object
spe <- readRDS(here::here("processed-data", "02_build_spe", "spe.rds"))

# Add variable of slide to colData(spe)
slide_df <- data.frame(spe$key, spe$sample_id)
slide_df <- slide_df %>%
    mutate(slide = case_when(
        grepl("Br1412", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br2706", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br3942", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br5242", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br5699", slide_df$spe.sample_id) ~ "V11D01387",
        grepl("Br6023", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br6129", slide_df$spe.sample_id) ~ "V11D01387",
        grepl("Br6299", slide_df$spe.sample_id) ~ "V11D01387",
        grepl("Br8195", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br8667", slide_df$spe.sample_id) ~ "V11L05334",
        grepl("Br8686", slide_df$spe.sample_id) ~ "V10B01088",
        grepl("Br8851", slide_df$spe.sample_id) ~ "V11D01387",
    ))

colData(spe)$slide <- factor(slide_df$slide, levels = c("V10B01088", "V11L05334", "V11D01387"))

## subset spe data based on slide
first_second_spe <- spe[, !spe$slide %in% c("V11D01387")]

vis_grid_clus(
    spe = first_second_spe,
    clustervar = "in_tissue",
    pdf = here::here("plots", "QC_plots", "First&Second_slides_in_tissue_grid.pdf"),
    sort_clust = FALSE,
    colors = c("TRUE" = "grey90", "FALSE" = "orange"),
    point_size = 2
)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(first_second_spe)$gene_name)
table(is_mito)
rowData(first_second_spe)$gene_name[is_mito]

# calculate per-spot QC metrics and store in colData
first_second_spe <- addPerCellQC(first_second_spe, subsets = list(mito = is_mito))
head(colData(first_second_spe))

# Create QC Metrics
qcstats <- perCellQCMetrics(first_second_spe, subsets = list(
    Mito = which(seqnames(first_second_spe) == "chrM")
))

qcfilter <- DataFrame(
    low_lib_size = isOutlier(qcstats$sum, type = "lower", log = TRUE, batch = first_second_spe$sample_id),
    low_n_features = isOutlier(qcstats$detected, type = "lower", log = TRUE, batch = first_second_spe$sample_id),
    high_subsets_Mito_percent = isOutlier(qcstats$subsets_Mito_percent, type = "higher", batch = first_second_spe$sample_id)
)

colSums(as.matrix(qcfilter))

qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent

# Add QC metrics to spe object
first_second_spe$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
first_second_spe$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
first_second_spe$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
first_second_spe$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

# plot histograms of QC metrics
pdf(file = here::here("plots", "QC_plots", "First&Second_slides_histograms_allSpots.pdf"), width = 8, height = 2.5)
par(mfrow = c(1, 3))
hist(colData(first_second_spe)$sum_umi, xlab = "sum", main = "UMIs per spot all donors")
hist(colData(first_second_spe)$sum_gene, xlab = "detected", main = "Genes per spot all donors")
hist(colData(first_second_spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito all donors")
par(mfrow = c(1, 1))
dev.off()

# QC plot of tissue spots discarded
for (i in colnames(qcfilter)) {
    vis_grid_clus(
        spe = first_second_spe,
        clustervar = paste0("scran_", i),
        pdf = here::here("plots", "QC_plots", paste0("First&Second_slides_scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange"),
        point_size = 2
    )
}

# remove combined set of low-quality spots
first_second_spe <- first_second_spe[, colData(first_second_spe)$scran_discard == FALSE]
dim(first_second_spe)

# calculate logcounts (log-transformed normalized counts) and store in object
first_second_spe$scran_quick_cluster <- quickCluster(first_second_spe, block = first_second_spe$sample_id)
first_second_spe <- computeSumFactors(first_second_spe, clusters = first_second_spe$scran_quick_cluster)
table(first_second_spe$scran_quick_cluster)
summary(sizeFactors(first_second_spe))
hist(sizeFactors(first_second_spe), breaks = 20)

first_second_spe <- logNormCounts(first_second_spe)

# check
assayNames(first_second_spe)

# Create directory for QC processed spe
dir_rdata <- here::here("processed-data", "QC_processed_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

saveRDS(first_second_spe, file = here::here("processed-data", "QC_processed_spe", "QCed_firsttwo_slides_spe.rds"))


## Reproducibility information
print("Reproducibility information:QC_spatial_DG_lifespan")
Sys.time()
proc.time()
options(width = 120)
session_info()
