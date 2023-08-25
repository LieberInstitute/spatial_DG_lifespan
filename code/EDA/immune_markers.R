########################################
# spatial_DG_lifespan project
# Plotting marker genes for immune cells
# Anthony Ramnauth, Aug 18 2023
########################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(ggnewscale)
    library(spatialLIBD)
    library(scater)
    library(scran)
    library(dplyr)
    library(sessioninfo)
})

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

rownames(spe) <- rowData(spe)$gene_name

## subset spe data based on BayesSpace clusters for DG
spe <- spe[, spe$bayesSpace_harmony_10 %in% c("1", "2", "4", "6", "7")]

bayes_df <- data.frame(spe$bayesSpace_harmony_10)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("1", spe.bayesSpace_harmony_10) ~ "SLM",
        grepl("2", spe.bayesSpace_harmony_10) ~ "ML",
        grepl("4", spe.bayesSpace_harmony_10) ~ "CA3&4",
        grepl("6", spe.bayesSpace_harmony_10) ~ "SGZ",
        grepl("7", spe.bayesSpace_harmony_10) ~ "GCL"
    ))

colData(spe)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("SLM", "ML", "CA3&4", "SGZ", "GCL"))

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

bay_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16",
    "SGZ" = "#1CFFCE", "GCL" = "#B00068")

features <- c(
    "C1QB",
    "CD83",
    "TREM2",
    "ID3",
    "P2RY12",
    "F13A1",
    "COLEC12",
    "LSP1",
    "SKAP1",
    "CD247",
    "CD163",
    "CD68",
    "MARCO",
    "MRC1",
    "MSR1",
    "FCGR3A",
    "KIR2DL4",
    "CD3E",
    "CD4",
    "CD8A",
    "FOXP3",
    "IL17A",
    "AIF1",
    "ITGAM"
)

pdf(file = here::here("plots","MHCII_genes", "violinplot_immune_genemarkers.pdf"))

plotExpression(spe, x = "BayesSpace", features = features, colour_by = "BayesSpace") +
    scale_color_manual(values  = bay_colors) +
    theme(plot.title = element_text(face = "italic"))

dev.off()

