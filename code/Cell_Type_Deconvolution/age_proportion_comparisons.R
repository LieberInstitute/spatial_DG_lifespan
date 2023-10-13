#####################################################
# spatial_DG_lifespan project
# Plotting Cell mean abundance changes with age
# Anthony Ramnauth, May 02 2023
#####################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(SpatialExperiment)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(spatialLIBD)
    library(viridis)
    library(ggplot2)
    library(ggcorrplot)
    library(ggsignif)
    library(ggh4x)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

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

# Isolate the entire DG to plot where cell types are most assigned
spe_DG <- spe[, spe$bayesSpace_harmony_10 == "2" |
        spe$bayesSpace_harmony_10 == "4" |
        spe$bayesSpace_harmony_10 == "6" |
        spe$bayesSpace_harmony_10 == "7"]

######################################
# Let's start with the DG locations
######################################

Bayes_df <-
    data.frame(spe_DG$key, spe_DG$sample_id, spe_DG$bayesSpace_harmony_10)
Bayes_df <- Bayes_df %>%
    mutate(
        BayesSpace = case_when(
            spe_DG.bayesSpace_harmony_10 == 2 ~ "ML",
            spe_DG.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe_DG.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe_DG.bayesSpace_harmony_10 == 7 ~ "GCL",
        )
    )

colData(spe_DG)$BayesSpace <-
    factor(Bayes_df$BayesSpace, levels = c("ML", "CA3&4", "SGZ", "GCL"))

# Create datafame of cell proportions and DG layer
cell_DGdf <- as.data.frame(colData(spe_DG)[, c(44:70)],
    row.names = spe_DG$key)

# isolate the abundance counts to convert to proportions

abundance_counts <- as.matrix(cell_DGdf[, c(1:25)])

# Convert mean abundances to proportions

convert_to_proportions <- function(row) {
  proportions <- row / sum(row)
  return(proportions)
}

# Apply the function to each row of the data frame
proportion_counts <- apply(abundance_counts, 1, convert_to_proportions)

proportion_counts <- as.data.frame(t(proportion_counts))
for (col in 1:ncol(proportion_counts)){
    colnames(proportion_counts)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(proportion_counts)[col])
}

stopifnot(rownames(proportion_counts) == rownames(cell_DGdf))

proportion_counts$age_bin <- cell_DGdf$age_bin
proportion_counts$BayesSpace <- cell_DGdf$BayesSpace

# Plot cell types of interest based on gene markers from DE results

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "DGlayer_age_proportions_cell2loc.pdf"), width = 16, height = 4)

strip <- strip_themed(background_x = elem_list_rect(fill = c("#E4E1E3", "#FEAF16", "#1CFFCE", "#B00068")))

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$Astro_1)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.45) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.45) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.45) +
    labs(x = "age", y = "Astro_1 estimated proportions") +
    theme_classic() +
    ylim(0, 0.5) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$Astro_2)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.75) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.75) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.75) +
    labs(x = "age", y = "Astro_2 estimated proportions") +
    theme_classic() +
    ylim(0, 0.8) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$InN_LAMP5)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.07) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.07) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.07) +
    labs(x = "age", y = "InN_LAMP5 estimated proportions") +
    theme_classic() +
    ylim(0, 0.1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$Macro)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.08) +
    labs(x = "age", y = "Macrophage estimated proportions") +
    theme_classic() +
    ylim(0, 0.1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$Microglia)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.08) +
    labs(x = "age", y = "Microglia estimated proportions") +
    theme_classic() +
    ylim(0, 0.1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$Myeloid)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.08) +
    labs(x = "age", y = "Myeloid estimated proportions") +
    theme_classic() +
    ylim(0, 0.1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$Oligo)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.4) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.4) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.4) +
    labs(x = "age", y = "Oligodendrocyte estimated proportions") +
    theme_classic() +
    ylim(0, 0.5) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(proportion_counts, aes(x = proportion_counts$age_bin, y = proportion_counts$T_Cell)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.08) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 0.08) +
    labs(x = "age", y = "T_Cell estimated proportions") +
    theme_classic() +
    ylim(0, 0.1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

dev.off()
