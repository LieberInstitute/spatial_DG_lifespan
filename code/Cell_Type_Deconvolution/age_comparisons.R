#####################################################
# spatial_DG_lifespan project
# Plotting Cell proportion changes with age with CARD
# Anthony Ramnauth, April 26 2023
#####################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(SpatialExperiment)
    library(dplyr)
    library(spatialLIBD)
    library(CARD)
    library(viridis)
    library(ggplot2)
    library(ggcorrplot)
    library(ggsignif)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_CARD_sestan.rds"))

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

# Separate out the DG layers of interest
spe_ML <- spe[, which(spe$bayesSpace_harmony_10 == "2")]
dim(spe_ML)

spe_CA3_4 <- spe[, which(spe$bayesSpace_harmony_10 == "4")]
dim(spe_CA3_4)

spe_SGZ <- spe[, which(spe$bayesSpace_harmony_10 == "6")]
dim(spe_SGZ)

spe_GCL <- spe[, which(spe$bayesSpace_harmony_10 == "7")]
dim(spe_GCL)

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
cell_DGdf <- as.data.frame(colData(spe_DG)[, c(44:73)],
    row.names = spe_DG$key)

cell_DGdf$dominant_cell_types <- NULL
cell_DGdf$age_bin <- NULL

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "DGlayer_cell.pdf"), width = 12, height = 8)

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$CA3_N)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "CA3 Proportion") +
    ggtitle("Cell-type Proportion of CA3") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Mossy)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Mossy Proportion") +
    ggtitle("Cell-type Proportion of Mossy") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$GC)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "GC Proportion") +
    ggtitle("Cell-type Proportion of GC") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_LAMP5)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_LAMP5 Proportion") +
    ggtitle("Cell-type Proportion of InN_LAMP5") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_PV)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_PV Proportion") +
    ggtitle("Cell-type Proportion of InN_PV") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_LHX6)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_LHX6 Proportion") +
    ggtitle("Cell-type Proportion of InN_LHX6") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_VIP)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_VIP Proportion") +
    ggtitle("Cell-type Proportion of InN_VIP") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_SST)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_SST Proportion") +
    ggtitle("Cell-type Proportion of InN_SST") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_NR2F2)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_NR2F2 Proportion") +
    ggtitle("Cell-type Proportion of InN_NR2F2") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_MEIS2)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_MEIS2 Proportion") +
    ggtitle("Cell-type Proportion of InN_MEIS2") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$InN_NPY)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "InN_NPY Proportion") +
    ggtitle("Cell-type Proportion of InN_NPY") +
    coord_cartesian(ylim=c(0, 0.002)) +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Oligo)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Oligo Proportion") +
    ggtitle("Cell-type Proportion of Oligo") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Microglia)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Microglia Proportion") +
    ggtitle("Cell-type Proportion of Microglia") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$OPC)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "OPC Proportion") +
    ggtitle("Cell-type Proportion of OPC") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Endoth)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Endoth Proportion") +
    ggtitle("Cell-type Proportion of Endoth") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$vSMC)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "vSMC Proportion") +
    ggtitle("Cell-type Proportion of vSMC") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$aSMC)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "aSMC Proportion") +
    ggtitle("Cell-type Proportion of aSMC") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$VLMC)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "VLMC Proportion") +
    ggtitle("Cell-type Proportion of VLMC") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Macro)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Macro Proportion") +
    ggtitle("Cell-type Proportion of Macro") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$COP)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "COP Proportion") +
    ggtitle("Cell-type Proportion of COP") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$T_Cell)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "T_Cell Proportion") +
    ggtitle("Cell-type Proportion of T_Cell") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Pericyte)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Pericyte Proportion") +
    ggtitle("Cell-type Proportion of Pericyte") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Myeloid)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Myeloid Proportion") +
    ggtitle("Cell-type Proportion of Myeloid") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Astro_1)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Astro_1 Proportion") +
    ggtitle("Cell-type Proportion of Astro_1") +
    theme_classic()

ggplot(cell_DGdf, aes(x = cell_DGdf$BayesSpace, y = cell_DGdf$Astro_2)) +
    geom_violin(aes(fill = factor(BayesSpace))) +
    geom_boxplot(width=0.1) +
    labs(x = "BayesSpace", y = "Astro_2 Proportion") +
    ggtitle("Cell-type Proportion of Astro_2") +
    theme_classic()

dev.off()


######################################
# Let's start with the Molecular Layer
######################################

# Create datafame of cell proportions and age_bin
cell_age_MLdf <- as.data.frame(colData(spe_ML)[, c(44:72)],
    age = spe_ML$age_bin,
    row.names = spe_ML$key)

cell_age_MLdf$dominant_cell_types <- NULL

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Age_changes_ML_cell.pdf"), width = 12, height = 8)

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_LAMP5)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_LAMP5 Proportion") +
    ggtitle("Cell-type Proportion of InN_LAMP5 in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_PV)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_PV Proportion") +
    ggtitle("Cell-type Proportion of InN_PV in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_LHX6)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_LHX6 Proportion") +
    ggtitle("Cell-type Proportion of InN_LHX6 in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_VIP)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_VIP Proportion") +
    ggtitle("Cell-type Proportion of InN_VIP in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_SST)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_SST Proportion") +
    ggtitle("Cell-type Proportion of InN_SST in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_NR2F2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    coord_cartesian(ylim=c(0, .06)) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.05) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.05) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.05) +
    labs(x = "age_bin", y = "InN_NR2F2 Proportion") +
    ggtitle("Cell-type Proportion of InN_NR2F2 in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_MEIS2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_MEIS2 Proportion") +
    ggtitle("Cell-type Proportion of InN_MEIS2 in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$InN_NPY)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_NPY Proportion") +
    ggtitle("Cell-type Proportion of InN_NPY in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Oligo)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Oligo Proportion") +
    ggtitle("Cell-type Proportion of Oligo in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Microglia)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Microglia Proportion") +
    ggtitle("Cell-type Proportion of Microglia in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$OPC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "OPC Proportion") +
    ggtitle("Cell-type Proportion of OPC in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Endoth)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Endoth Proportion") +
    ggtitle("Cell-type Proportion of Endoth in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$vSMC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "vSMC Proportion") +
    ggtitle("Cell-type Proportion of vSMC in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$VLMC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "VLMC Proportion") +
    ggtitle("Cell-type Proportion of VLMC in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Macro)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Macro Proportion") +
    ggtitle("Cell-type Proportion of Macro in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$COP)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "COP Proportion") +
    ggtitle("Cell-type Proportion of COP in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$T_Cell)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "T_Cell Proportion") +
    ggtitle("Cell-type Proportion of T_Cell in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Pericyte)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Pericyte Proportion") +
    ggtitle("Cell-type Proportion of Pericyte in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Myeloid)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Myeloid Proportion") +
    ggtitle("Cell-type Proportion of Myeloid in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Astro_1)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Astro_1 Proportion") +
    ggtitle("Cell-type Proportion of Astro_1 in ML") +
    theme_classic()

ggplot(cell_age_MLdf, aes(x = cell_age_MLdf$age_bin, y = cell_age_MLdf$Astro_2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Astro_2 Proportion") +
    ggtitle("Cell-type Proportion of Astro_2 in ML") +
    theme_classic()

dev.off()

######################
# Let's move on to GCL
######################

# Create datafame of cell proportions and age_bin
cell_age_GCLdf <- as.data.frame(colData(spe_GCL)[, c(44:72)],
    age = spe_GCL$age_bin,
    row.names = spe_GCL$key)

cell_age_GCLdf$dominant_cell_types <- NULL

# Loop through each numeric column variable and create a violin plot

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Age_changes_GCL_cell.pdf"), width = 12, height = 8)

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$GC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "GC Proportion") +
    ggtitle("Cell-type Proportion of GC in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_LAMP5)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_LAMP5 Proportion") +
    ggtitle("Cell-type Proportion of InN_LAMP5 in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_PV)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_PV Proportion") +
    ggtitle("Cell-type Proportion of InN_PV in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_LHX6)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_LHX6 Proportion") +
    ggtitle("Cell-type Proportion of InN_LHX6 in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_VIP)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_VIP Proportion") +
    ggtitle("Cell-type Proportion of InN_VIP in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_SST)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_SST Proportion") +
    ggtitle("Cell-type Proportion of InN_SST in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_NR2F2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    coord_cartesian(ylim=c(0, .06)) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.05) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.05) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.05) +
    labs(x = "age_bin", y = "InN_NR2F2 Proportion") +
    ggtitle("Cell-type Proportion of InN_NR2F2 in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_MEIS2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_MEIS2 Proportion") +
    ggtitle("Cell-type Proportion of InN_MEIS2 in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$InN_NPY)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_NPY Proportion") +
    ggtitle("Cell-type Proportion of InN_NPY in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Oligo)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Oligo Proportion") +
    ggtitle("Cell-type Proportion of Oligo in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Microglia)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Microglia Proportion") +
    ggtitle("Cell-type Proportion of Microglia in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$OPC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "OPC Proportion") +
    ggtitle("Cell-type Proportion of OPC in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Endoth)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Endoth Proportion") +
    ggtitle("Cell-type Proportion of Endoth in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$vSMC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "vSMC Proportion") +
    ggtitle("Cell-type Proportion of vSMC in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$VLMC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "VLMC Proportion") +
    ggtitle("Cell-type Proportion of VLMC in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Macro)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Macro Proportion") +
    ggtitle("Cell-type Proportion of Macro in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$COP)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "COP Proportion") +
    ggtitle("Cell-type Proportion of COP in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$T_Cell)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "T_Cell Proportion") +
    ggtitle("Cell-type Proportion of T_Cell in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Pericyte)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Pericyte Proportion") +
    ggtitle("Cell-type Proportion of Pericyte in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Myeloid)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Myeloid Proportion") +
    ggtitle("Cell-type Proportion of Myeloid in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Astro_1)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Astro_1 Proportion") +
    ggtitle("Cell-type Proportion of Astro_1 in GCL") +
    theme_classic()

ggplot(cell_age_GCLdf, aes(x = cell_age_GCLdf$age_bin, y = cell_age_GCLdf$Astro_2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Astro_2 Proportion") +
    ggtitle("Cell-type Proportion of Astro_2 in GCL") +
    theme_classic()

dev.off()

######################
# Let's move on to SGZ
######################

# Create datafame of cell proportions and age_bin
cell_age_SGZdf <- as.data.frame(colData(spe_SGZ)[, c(44:72)],
    age = spe_SGZ$age_bin,
    row.names = spe_SGZ$key)

cell_age_SGZdf$dominant_cell_types <- NULL

# Loop through each numeric column variable and create a violin plot

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "Age_changes_SGZ_cell.pdf"), width = 12, height = 8)

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_LAMP5)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_LAMP5 Proportion") +
    ggtitle("Cell-type Proportion of InN_LAMP5 in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_PV)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_PV Proportion") +
    ggtitle("Cell-type Proportion of InN_PV in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_LHX6)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_LHX6 Proportion") +
    ggtitle("Cell-type Proportion of InN_LHX6 in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_VIP)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_VIP Proportion") +
    ggtitle("Cell-type Proportion of InN_VIP in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_SST)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_SST Proportion") +
    ggtitle("Cell-type Proportion of InN_SST in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_NR2F2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1, outlier.shape = NA) +
    coord_cartesian(ylim=c(0, .06)) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE, y_position = 0.05) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE, y_position = 0.05) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE, y_position = 0.05) +
    labs(x = "age_bin", y = "InN_NR2F2 Proportion") +
    ggtitle("Cell-type Proportion of InN_NR2F2 in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_MEIS2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_MEIS2 Proportion") +
    ggtitle("Cell-type Proportion of InN_MEIS2 in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$InN_NPY)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "InN_NPY Proportion") +
    ggtitle("Cell-type Proportion of InN_NPY in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Oligo)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Oligo Proportion") +
    ggtitle("Cell-type Proportion of Oligo in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Microglia)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Microglia Proportion") +
    ggtitle("Cell-type Proportion of Microglia in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$OPC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "OPC Proportion") +
    ggtitle("Cell-type Proportion of OPC in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Endoth)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Endoth Proportion") +
    ggtitle("Cell-type Proportion of Endoth in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$vSMC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "vSMC Proportion") +
    ggtitle("Cell-type Proportion of vSMC in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$VLMC)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "VLMC Proportion") +
    ggtitle("Cell-type Proportion of VLMC in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Macro)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Macro Proportion") +
    ggtitle("Cell-type Proportion of Macro in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$COP)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "COP Proportion") +
    ggtitle("Cell-type Proportion of COP in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$T_Cell)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "T_Cell Proportion") +
    ggtitle("Cell-type Proportion of T_Cell in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Pericyte)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Pericyte Proportion") +
    ggtitle("Cell-type Proportion of Pericyte in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Myeloid)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Myeloid Proportion") +
    ggtitle("Cell-type Proportion of Myeloid in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Astro_1)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Astro_1 Proportion") +
    ggtitle("Cell-type Proportion of Astro_1 in SGZ") +
    theme_classic()

ggplot(cell_age_SGZdf, aes(x = cell_age_SGZdf$age_bin, y = cell_age_SGZdf$Astro_2)) +
    geom_violin(aes(fill = factor(age_bin))) +
    geom_boxplot(width=0.1) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Infant", "Teen")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Teen", "Adult")),
              map_signif_level=TRUE) +
    geom_signif(test = "wilcox.test", comparisons = list(c("Adult", "Elderly")),
              map_signif_level=TRUE) +
    labs(x = "age_bin", y = "Astro_2 Proportion") +
    ggtitle("Cell-type Proportion of Astro_2 in SGZ") +
    theme_classic()

dev.off()


sessionInfo()
