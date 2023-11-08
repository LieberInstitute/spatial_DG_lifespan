###############################
# spatial_DG_lifespan project
# Plotting cell counts per spot
# Anthony Ramnauth, Nov 06 2023
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(ggh4x)
    library(dplyr)
})

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Retain the capture areas with high quality H&E images to assess cell segmentation

spe_subset <- spe[, spe$sample_id == "Br8686" |
        spe$sample_id == "Br2706" |
        spe$sample_id == "Br3942" |
        spe$sample_id == "Br6023" |
        spe$sample_id == "Br8195" |
        spe$sample_id == "Br1412" |
        spe$sample_id == "Br8667" |
        spe$sample_id == "Br5242"]

# Setting up data for  plots of counts vs age faceted by spatial domain

df <-
    data.frame(spe_subset$key, spe_subset$sample_id, spe_subset$bayesSpace_harmony_10)

df <- df %>%
    mutate(
        BayesSpace = case_when(
            spe_subset.bayesSpace_harmony_10 == 1 ~ "SLM",
            spe_subset.bayesSpace_harmony_10 == 2 ~ "ML",
            spe_subset.bayesSpace_harmony_10 == 3 ~ "CP",
            spe_subset.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe_subset.bayesSpace_harmony_10 == 5 ~ "SR",
            spe_subset.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe_subset.bayesSpace_harmony_10 == 7 ~ "GCL",
            spe_subset.bayesSpace_harmony_10 == 8 ~ "SL",
            spe_subset.bayesSpace_harmony_10 == 9 ~ "CA1",
            spe_subset.bayesSpace_harmony_10 == 10 ~ "WM",
        )
    )

colData(spe_subset)$BayesSpace <-
    factor(df$BayesSpace, levels = c("WM", "SLM", "SL", "ML", "CP", "SR", "SGZ",
        "CA3&4", "CA1", "GCL"))

df <- as.data.frame(colData(spe_subset))

# Remove CP cluster since that has higher variance and was not included in upstream DE analyses
# Remove Choroid Plexus cluster
df$spatialnum <- spe_subset$bayesSpace_harmony_10
df = df[which(df$spatialnum != "3"), ]

strip <- strip_themed(background_x = elem_list_rect(fill = c("#2ED9FF", "#5A5156", "#90AD1C", "#E4E1E3",
    "#FE00FA", "#1CFFCE", "#FEAF16", "#16FF32", "#B00068")))

pdf(file = here::here("plots", "Trajectories", "CellperSpot_vs_age.pdf"), width = 5, height = 5)

ggplot(df, aes(x = age, y = NBW)) +
    geom_boxplot(width = 10, notch = F, outlier.colour = NA,  mapping = aes(fill = as.factor(age))) +
    geom_smooth(data = df, mapping = aes(x = age, y = NBW), color='red3', se = F, method = 'loess') +
    labs(x = "age", y = "Est. Cell count") +
    facet_wrap2(~ BayesSpace, nrow = 4, ncol = 3, strip = strip) +
    theme(axis.text.x = element_text(size = 7), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank()) +
    guides(fill = FALSE) +
    theme_classic()

ggplot(df, aes(x = age, y = PBW)) +
    geom_boxplot(width = 10, notch = F, outlier.colour = NA,  mapping = aes(fill = as.factor(age))) +
    geom_smooth(data = df, mapping = aes(x = age, y = PBW), color='red3', se = F, method = 'loess') +
    labs(x = "age", y = "% spot covered with cells") +
    facet_wrap2(~ BayesSpace, nrow = 4, ncol = 3, strip = strip) +
    theme(axis.text.x = element_text(size = 7), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank()) +
    guides(fill = FALSE) +
    theme_classic()

dev.off()

