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
    library(spatialLIBD)
    library(viridis)
    library(ggplot2)
    library(ggcorrplot)
    library(ggsignif)
    library(ggh4x)
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

# Isolate SLM (absent cluster 10 WM) for cell type abundance

spe_SLM <- spe[, spe$bayesSpace_harmony_10 == "1"]

spe_SLM$bayesSpace_harmony_10 <- as.numeric(spe_SLM$bayesSpace_harmony_10)
spe_SLM$bayesSpace_harmony_10 <- as.factor(spe_SLM$bayesSpace_harmony_10)

SLM_df <-
    data.frame(spe_SLM$key, spe_SLM$sample_id, spe_SLM$bayesSpace_harmony_10)
SLM_df <- SLM_df %>%
    mutate(
        BayesSpace = case_when(
            spe_SLM.bayesSpace_harmony_10 == 1 ~ "SLM"
        )
    )

colData(spe_SLM)$BayesSpace <-
    factor(SLM_df$BayesSpace, levels = c("SLM"))

# Create datafame of cell proportions and DG layer
cell_SLMdf <- as.data.frame(colData(spe_SLM)[, c(44:70)],
    row.names = spe_SLM$key)

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


##################################
# Facet to see everything at once.
##################################

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "DGlayer_age_facet_cell2loc.pdf"), width = 16, height = 4)

strip <- strip_themed(background_x = elem_list_rect(fill = c("#E4E1E3", "#FEAF16", "#1CFFCE", "#B00068")))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Astro_1)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 4) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 4) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE,  y_position = 4) +
    labs(x = "age", y = "Astro_1 mean abundance") +
    ggtitle("Mean abundance of Astro_1") +
    theme_classic() +
    ylim(0, 5) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Astro_2)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 6) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 6) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 6) +
    labs(x = "age", y = "Astro_2 mean abundance") +
    ggtitle("Mean abundance of Astro_2") +
    theme_classic() +
    ylim(0, 7) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_CA1_N)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "CA1_N mean abundance") +
    ggtitle("Mean abundance of CA1_N") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_CA2_N)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "CA2_N mean abundance") +
    ggtitle("Mean abundance of CA2_N") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_CA3_N)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "CA3_N mean abundance") +
    ggtitle("Mean abundance of CA3_N") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_COP)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "COP mean abundance") +
    ggtitle("Mean abundance of COP") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Endoth)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Endoth mean abundance") +
    ggtitle("Mean abundance of Endoth") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_GC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "GC mean abundance") +
    ggtitle("Mean abundance of GC") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_LAMP5)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.5) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.5) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.5) +
    labs(x = "age", y = "InN_LAMP5 mean abundance") +
    ggtitle("Mean abundance of InN_LAMP5") +
    theme_classic() +
    ylim(0, 1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_LHX6)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_LHX6 mean abundance") +
    ggtitle("Mean abundance of InN_LHX6") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))


ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_MEIS2)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_MEIS2 mean abundance") +
    ggtitle("Mean abundance of InN_MEIS2") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_NR2F2)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_NR2F2 mean abundance") +
    ggtitle("Mean abundance of InN_NR2F2") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_PV)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_PV mean abundance") +
    ggtitle("Mean abundance of InN_PV") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_SST)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 4) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 4) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 4) +
    labs(x = "age", y = "InN_SST mean abundance") +
    ggtitle("Mean abundance of InN_SST") +
    theme_classic() +
    ylim(0, 5) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_InN_VIP)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_VIP mean abundance") +
    ggtitle("Mean abundance of InN_VIP") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Macro)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.8) +
    labs(x = "age", y = "Macrophage mean abundance") +
    ggtitle("Mean abundance of Macrophages") +
    theme_classic() +
    ylim(0, 1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Microglia)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.8) +
    labs(x = "age", y = "Microglia mean abundance") +
    ggtitle("Mean abundance of Microglia") +
    theme_classic() +
    ylim(0, 1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Mossy)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Mossy mean abundance") +
    ggtitle("Mean abundance of Mossy") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Myeloid)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.8) +
    labs(x = "age", y = "Myeloid mean abundance") +
    ggtitle("Mean abundance of Myeloid cells") +
    theme_classic() +
    ylim(0, 1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_OPC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "OPC mean abundance") +
    ggtitle("Mean abundance of OPC") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Oligo)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 3) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 3) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 3) +
    labs(x = "age", y = "Oligodendrocyte mean abundance") +
    ggtitle("Mean abundance of Oligodendrocytes") +
    theme_classic() +
    ylim(0, 4) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_Pericyte)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Pericyte mean abundance") +
    ggtitle("Mean abundance of Pericyte") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_SMC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "SMC mean abundance") +
    ggtitle("Mean abundance of SMC") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_T_Cell)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.8) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.8) +
    labs(x = "age", y = "T cell mean abundance") +
    ggtitle("Mean abundance of T cells") +
    theme_classic() +
    ylim(0, 1) +
    facet_wrap2(~ BayesSpace, nrow = 1, ncol = 4, strip = strip) +
    theme(text = element_text(size = 16), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(cell_DGdf, aes(x = cell_DGdf$age_bin, y = cell_DGdf$meanscell_abundance_w_sf_VLMC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "VLMC mean abundance") +
    ggtitle("Mean abundance of VLMC") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

dev.off()

##################################################################################################

# Plot for SLM only

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "SLM_age_cell2loc.pdf"), width = 8, height = 8)

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Astro_1)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 5) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 5) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 5) +
    labs(x = "age", y = "Astro_1 mean abundance") +
    ggtitle("Mean abundance of Astro_1") +
    theme_classic() +
    ylim(0, 6) +
    theme(text = element_text(size = 20))

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Astro_2)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 7) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 7) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 7) +
    labs(x = "age", y = "Astro_2 mean abundance") +
    ggtitle("Mean abundance of Astro_2") +
    theme_classic() +
    ylim(0, 8) +
    theme(text = element_text(size = 20))

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_CA1_N)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "CA1_N mean abundance") +
    ggtitle("Mean abundance of CA1_N") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_CA2_N)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "CA2_N mean abundance") +
    ggtitle("Mean abundance of CA2_N") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_CA3_N)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "CA3_N mean abundance") +
    ggtitle("Mean abundance of CA3_N") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_COP)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "COP mean abundance") +
    ggtitle("Mean abundance of COP") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Endoth)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Endoth mean abundance") +
    ggtitle("Mean abundance of Endoth") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_GC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "GC mean abundance") +
    ggtitle("Mean abundance of GC") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_LAMP5)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_LAMP5 mean abundance") +
    ggtitle("Mean abundance of InN_LAMP5") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_LHX6)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_LHX6 mean abundance") +
    ggtitle("Mean abundance of InN_LHX6") +
    theme_classic()


ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_MEIS2)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_MEIS2 mean abundance") +
    ggtitle("Mean abundance of InN_MEIS2") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_NR2F2)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_NR2F2 mean abundance") +
    ggtitle("Mean abundance of InN_NR2F2") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_PV)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_PV mean abundance") +
    ggtitle("Mean abundance of InN_PV") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_SST)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_SST mean abundance") +
    ggtitle("Mean abundance of InN_SST") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_InN_VIP)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "InN_VIP mean abundance") +
    ggtitle("Mean abundance of InN_VIP") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Macro)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Macro mean abundance") +
    ggtitle("Mean abundance of Macro") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Microglia)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.85) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.85) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.85) +
    labs(x = "age", y = "Microglia mean abundance") +
    ggtitle("Mean abundance of Microglia") +
    theme_classic() +
    ylim(0, 1) +
    theme(text = element_text(size = 20))

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Mossy)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Mossy mean abundance") +
    ggtitle("Mean abundance of Mossy") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Myeloid)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Myeloid mean abundance") +
    ggtitle("Mean abundance of Myeloid") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_OPC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "OPC mean abundance") +
    ggtitle("Mean abundance of OPC") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Oligo)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Oligo mean abundance") +
    ggtitle("Mean abundance of Oligo") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_Pericyte)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "Pericyte mean abundance") +
    ggtitle("Mean abundance of Pericyte") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_SMC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "SMC mean abundance") +
    ggtitle("Mean abundance of SMC") +
    theme_classic() +
    facet_wrap(vars(BayesSpace))

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_T_Cell)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "T_Cell mean abundance") +
    ggtitle("Mean abundance of T_Cell") +
    theme_classic()

ggplot(cell_SLMdf, aes(x = cell_SLMdf$age_bin, y = cell_SLMdf$meanscell_abundance_w_sf_VLMC)) +
    geom_violin(aes(fill = age_bin)) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age", y = "VLMC mean abundance") +
    ggtitle("Mean abundance of VLMC") +
    theme_classic()

dev.off()
