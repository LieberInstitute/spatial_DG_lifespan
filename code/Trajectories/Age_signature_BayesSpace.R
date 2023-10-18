###################################################################
# spatial_DG_lifespan project
# Common Age signature analysis controlling for BayesSpace clusters
# Anthony Ramnauth, Oct 16 2023
###################################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(sessioninfo)
    library(SingleCellExperiment)
    library(rafalib)
    library(limma)
    library(edgeR)
    library(scran)
    library(VISION)
    library(dplyr)
    library(ggh4x)
    library(ggsignif)
    library(spatialLIBD)
    library(lsmeans)
    library(RColorBrewer)
})

# Load .csv files from DE of infant vs. age_bins

Teen1 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace1_DE.csv"))
Teen1up <- Teen1[Teen1$logFC >= 1.5,]
Teen1down <- Teen1[Teen1$logFC <= -1.5,]

Teen2 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace2_DE.csv"))
Teen2up <- Teen2[Teen2$logFC >= 1.5,]
Teen2down <- Teen2[Teen2$logFC <= -1.5,]

Teen4 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace4_DE.csv"))
Teen4up <- Teen4[Teen4$logFC >= 1.5,]
Teen4down <- Teen4[Teen4$logFC <= -1.5,]

Teen5 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace5_DE.csv"))
Teen5up <- Teen5[Teen5$logFC >= 1.5,]
Teen5down <- Teen5[Teen5$logFC <= -1.5,]

Teen6 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace6_DE.csv"))
Teen6up <- Teen6[Teen6$logFC >= 1.5,]
Teen6down <- Teen6[Teen6$logFC <= -1.5,]

Teen7 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace7_DE.csv"))
Teen7up <- Teen7[Teen7$logFC >= 1.5,]
Teen7down <- Teen7[Teen7$logFC <= -1.5,]

Teen8 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace8_DE.csv"))
Teen8up <- Teen8[Teen8$logFC >= 1.5,]
Teen8down <- Teen8[Teen8$logFC <= -1.5,]

Teen9 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace9_DE.csv"))
Teen9up <- Teen9[Teen9$logFC >= 1.5,]
Teen9down <- Teen9[Teen9$logFC <= -1.5,]

Teen10 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "TeenvsInfant_BayesSpace10_DE.csv"))
Teen10up <- Teen10[Teen10$logFC >= 1.5,]
Teen10down <- Teen10[Teen10$logFC <= -1.5,]

Adult1 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace1_DE.csv"))
Adult1up <- Adult1[Adult1$logFC >= 1.5,]
Adult1down <- Adult1[Adult1$logFC <= -1.5,]

Adult2 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace2_DE.csv"))
Adult2up <- Adult2[Adult2$logFC >= 1.5,]
Adult2down <- Adult2[Adult2$logFC <= -1.5,]

Adult4 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace4_DE.csv"))
Adult4up <- Adult4[Adult4$logFC >= 1.5,]
Adult4down <- Adult4[Adult4$logFC <= -1.5,]

Adult5 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace5_DE.csv"))
Adult5up <- Adult5[Adult5$logFC >= 1.5,]
Adult5down <- Adult5[Adult5$logFC <= -1.5,]

Adult6 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace6_DE.csv"))
Adult6up <- Adult6[Adult6$logFC >= 1.5,]
Adult6down <- Adult6[Adult6$logFC <= -1.5,]

Adult7 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace7_DE.csv"))
Adult7up <- Adult7[Adult7$logFC >= 1.5,]
Adult7down <- Adult7[Adult7$logFC <= -1.5,]

Adult8 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace8_DE.csv"))
Adult8up <- Adult8[Adult8$logFC >= 1.5,]
Adult8down <- Adult8[Adult8$logFC <= -1.5,]

Adult9 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace9_DE.csv"))
Adult9up <- Adult9[Adult9$logFC >= 1.5,]
Adult9down <- Adult9[Adult9$logFC <= -1.5,]

Adult10 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "AdultvsInfant_BayesSpace10_DE.csv"))
Adult10up <- Adult10[Adult10$logFC >= 1.5,]
Adult10down <- Adult10[Adult10$logFC <= -1.5,]

Elderly1 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace1_DE.csv"))
Elderly1up <- Elderly1[Elderly1$logFC >= 1.5,]
Elderly1down <- Elderly1[Elderly1$logFC <= -1.5,]

Elderly2 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace2_DE.csv"))
Elderly2up <- Elderly2[Elderly2$logFC >= 1.5,]
Elderly2down <- Elderly2[Elderly2$logFC <= -1.5,]

Elderly4 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace4_DE.csv"))
Elderly4up <- Elderly4[Elderly4$logFC >= 1.5,]
Elderly4down <- Elderly4[Elderly4$logFC <= -1.5,]

Elderly5 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace5_DE.csv"))
Elderly5up <- Elderly5[Elderly5$logFC >= 1.5,]
Elderly5down <- Elderly5[Elderly5$logFC <= -1.5,]

Elderly6 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace6_DE.csv"))
Elderly6up <- Elderly6[Elderly6$logFC >= 1.5,]
Elderly6down <- Elderly6[Elderly6$logFC <= -1.5,]

Elderly7 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace7_DE.csv"))
Elderly7up <- Elderly7[Elderly7$logFC >= 1.5,]
Elderly7down <- Elderly7[Elderly7$logFC <= -1.5,]

Elderly8 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace8_DE.csv"))
Elderly8up <- Elderly8[Elderly8$logFC >= 1.5,]
Elderly8down <- Elderly8[Elderly8$logFC <= -1.5,]

Elderly9 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace9_DE.csv"))
Elderly9up <- Elderly9[Elderly9$logFC >= 1.5,]
Elderly9down <- Elderly9[Elderly9$logFC <= -1.5,]

Elderly10 <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results",
    "ElderlyvsInfant_BayesSpace10_DE.csv"))
Elderly10up <- Elderly10[Elderly10$logFC >= 1.5,]
Elderly10down <- Elderly10[Elderly10$logFC <= -1.5,]

# Check, within each BayesSpace cluster, whether gene is significant for >= 2 DE tests

#######################
# Increasing expression
#######################

# SLM
combo1_up_SLM <- Teen1up[Teen1up$gene_id %in% Adult1up$gene_id,]
combo2_up_SLM <- Teen1up[Teen1up$gene_id %in% Elderly1up$gene_id,]
combo3_up_SLM <- Adult1up[Adult1up$gene_id %in% Elderly1up$gene_id,]

# Combine and find unique genes
combo_up_SLM <- rbind(combo1_up_SLM, combo2_up_SLM, combo3_up_SLM)
combo_up_SLM <- unique(combo_up_SLM[,1:3])

# ML
combo1_up_ML <- Teen2up[Teen2up$gene_id %in% Adult2up$gene_id,]
combo2_up_ML <- Teen2up[Teen2up$gene_id %in% Elderly2up$gene_id,]
combo3_up_ML <- Adult2up[Adult2up$gene_id %in% Elderly2up$gene_id,]

# Combine and find unique genes
combo_up_ML <- rbind(combo1_up_ML, combo2_up_ML, combo3_up_ML)
combo_up_ML <- unique(combo_up_ML[,1:3])

# CA3&4
combo1_up_CA3_4 <- Teen4up[Teen4up$gene_id %in% Adult4up$gene_id,]
combo2_up_CA3_4 <- Teen4up[Teen4up$gene_id %in% Elderly4up$gene_id,]
combo3_up_CA3_4 <- Adult4up[Adult4up$gene_id %in% Elderly4up$gene_id,]

# Combine and find unique genes
combo_up_CA3_4 <- rbind(combo1_up_CA3_4, combo2_up_CA3_4, combo3_up_CA3_4)
combo_up_CA3_4 <- unique(combo_up_CA3_4[,1:3])

# SR
combo1_up_SR <- Teen5up[Teen5up$gene_id %in% Adult5up$gene_id,]
combo2_up_SR <- Teen5up[Teen5up$gene_id %in% Elderly5up$gene_id,]
combo3_up_SR <- Adult5up[Adult5up$gene_id %in% Elderly5up$gene_id,]

# Combine and find unique genes
combo_up_SR <- rbind(combo1_up_SR, combo2_up_SR, combo3_up_SR)
combo_up_SR <- unique(combo_up_SR[,1:3])

# SGZ
combo1_up_SGZ <- Teen6up[Teen6up$gene_id %in% Adult6up$gene_id,]
combo2_up_SGZ <- Teen6up[Teen6up$gene_id %in% Elderly6up$gene_id,]
combo3_up_SGZ <- Adult6up[Adult6up$gene_id %in% Elderly6up$gene_id,]

# Combine and find unique genes
combo_up_SGZ <- rbind(combo1_up_SGZ, combo2_up_SGZ, combo3_up_SGZ)
combo_up_SGZ <- unique(combo_up_SGZ[,1:3])

# GCL
combo1_up_GCL <- Teen7up[Teen7up$gene_id %in% Adult7up$gene_id,]
combo2_up_GCL <- Teen7up[Teen7up$gene_id %in% Elderly7up$gene_id,]
combo3_up_GCL <- Adult7up[Adult7up$gene_id %in% Elderly7up$gene_id,]

# Combine and find unique genes
combo_up_GCL <- rbind(combo1_up_GCL, combo2_up_GCL, combo3_up_GCL)
combo_up_GCL <- unique(combo_up_GCL[,1:3])

# SL
combo1_up_SL <- Teen8up[Teen8up$gene_id %in% Adult8up$gene_id,]
combo2_up_SL <- Teen8up[Teen8up$gene_id %in% Elderly8up$gene_id,]
combo3_up_SL <- Adult8up[Adult8up$gene_id %in% Elderly8up$gene_id,]

# Combine and find unique genes
combo_up_SL <- rbind(combo1_up_SL, combo2_up_SL, combo3_up_SL)
combo_up_SL <- unique(combo_up_SL[,1:3])

# CA1
combo1_up_CA1 <- Teen9up[Teen9up$gene_id %in% Adult9up$gene_id,]
combo2_up_CA1 <- Teen9up[Teen9up$gene_id %in% Elderly9up$gene_id,]
combo3_up_CA1 <- Adult9up[Adult9up$gene_id %in% Elderly9up$gene_id,]

# Combine and find unique genes
combo_up_CA1 <- rbind(combo1_up_CA1, combo2_up_CA1, combo3_up_CA1)
combo_up_CA1 <- unique(combo_up_CA1[,1:3])

# WM
combo1_up_WM <- Teen10up[Teen10up$gene_id %in% Adult10up$gene_id,]
combo2_up_WM <- Teen10up[Teen10up$gene_id %in% Elderly10up$gene_id,]
combo3_up_WM <- Adult10up[Adult10up$gene_id %in% Elderly10up$gene_id,]

# Combine and find unique genes
combo_up_WM <- rbind(combo1_up_WM, combo2_up_WM, combo3_up_WM)
combo_up_WM <- unique(combo_up_WM[,1:3])

Increasing <- rbind(combo_up_CA1, combo_up_CA3_4, combo_up_GCL, combo_up_ML, combo_up_SGZ,
    combo_up_SL, combo_up_SLM, combo_up_SR, combo_up_WM)

# Isolate rows repeated > 2 times for stringency in positive aging signature

Common_increasing <- Increasing %>%
  group_by_all() %>%
  filter(n() > 2) %>%
  as.data.frame() %>%
  unique()

#######################
# Decreasing expression
#######################

# SLM
combo1_down_SLM <- Teen1down[Teen1down$gene_id %in% Adult1down$gene_id,]
combo2_down_SLM <- Teen1down[Teen1down$gene_id %in% Elderly1down$gene_id,]
combo3_down_SLM <- Adult1down[Adult1down$gene_id %in% Elderly1down$gene_id,]

# Combine and find unique genes
combo_down_SLM <- rbind(combo1_down_SLM, combo2_down_SLM, combo3_down_SLM)
combo_down_SLM <- unique(combo_down_SLM[,1:3])

# ML
combo1_down_ML <- Teen2down[Teen2down$gene_id %in% Adult2down$gene_id,]
combo2_down_ML <- Teen2down[Teen2down$gene_id %in% Elderly2down$gene_id,]
combo3_down_ML <- Adult2down[Adult2down$gene_id %in% Elderly2down$gene_id,]

# Combine and find unique genes
combo_down_ML <- rbind(combo1_down_ML, combo2_down_ML, combo3_down_ML)
combo_down_ML <- unique(combo_down_ML[,1:3])

# CA3&4
combo1_down_CA3_4 <- Teen4down[Teen4down$gene_id %in% Adult4down$gene_id,]
combo2_down_CA3_4 <- Teen4down[Teen4down$gene_id %in% Elderly4down$gene_id,]
combo3_down_CA3_4 <- Adult4down[Adult4down$gene_id %in% Elderly4down$gene_id,]

# Combine and find unique genes
combo_down_CA3_4 <- rbind(combo1_down_CA3_4, combo2_down_CA3_4, combo3_down_CA3_4)
combo_down_CA3_4 <- unique(combo_down_CA3_4[,1:3])

# SR
combo1_down_SR <- Teen5down[Teen5down$gene_id %in% Adult5down$gene_id,]
combo2_down_SR <- Teen5down[Teen5down$gene_id %in% Elderly5down$gene_id,]
combo3_down_SR <- Adult5down[Adult5down$gene_id %in% Elderly5down$gene_id,]

# Combine and find unique genes
combo_down_SR <- rbind(combo1_down_SR, combo2_down_SR, combo3_down_SR)
combo_down_SR <- unique(combo_down_SR[,1:3])

# SGZ
combo1_down_SGZ <- Teen6down[Teen6down$gene_id %in% Adult6down$gene_id,]
combo2_down_SGZ <- Teen6down[Teen6down$gene_id %in% Elderly6down$gene_id,]
combo3_down_SGZ <- Adult6down[Adult6down$gene_id %in% Elderly6down$gene_id,]

# Combine and find unique genes
combo_down_SGZ <- rbind(combo1_down_SGZ, combo2_down_SGZ, combo3_down_SGZ)
combo_down_SGZ <- unique(combo_down_SGZ[,1:3])

# GCL
combo1_down_GCL <- Teen7down[Teen7down$gene_id %in% Adult7down$gene_id,]
combo2_down_GCL <- Teen7down[Teen7down$gene_id %in% Elderly7down$gene_id,]
combo3_down_GCL <- Adult7down[Adult7down$gene_id %in% Elderly7down$gene_id,]

# Combine and find unique genes
combo_down_GCL <- rbind(combo1_down_GCL, combo2_down_GCL, combo3_down_GCL)
combo_down_GCL <- unique(combo_down_GCL[,1:3])

# SL
combo1_down_SL <- Teen8down[Teen8down$gene_id %in% Adult8down$gene_id,]
combo2_down_SL <- Teen8down[Teen8down$gene_id %in% Elderly8down$gene_id,]
combo3_down_SL <- Adult8down[Adult8down$gene_id %in% Elderly8down$gene_id,]

# Combine and find unique genes
combo_down_SL <- rbind(combo1_down_SL, combo2_down_SL, combo3_down_SL)
combo_down_SL <- unique(combo_down_SL[,1:3])

# CA1
combo1_down_CA1 <- Teen9down[Teen9down$gene_id %in% Adult9down$gene_id,]
combo2_down_CA1 <- Teen9down[Teen9down$gene_id %in% Elderly9down$gene_id,]
combo3_down_CA1 <- Adult9down[Adult9down$gene_id %in% Elderly9down$gene_id,]

# Combine and find unique genes
combo_down_CA1 <- rbind(combo1_down_CA1, combo2_down_CA1, combo3_down_CA1)
combo_down_CA1 <- unique(combo_down_CA1[,1:3])

# WM
combo1_down_WM <- Teen10down[Teen10down$gene_id %in% Adult10down$gene_id,]
combo2_down_WM <- Teen10down[Teen10down$gene_id %in% Elderly10down$gene_id,]
combo3_down_WM <- Adult10down[Adult10down$gene_id %in% Elderly10down$gene_id,]

# Combine and find unique genes
combo_down_WM <- rbind(combo1_down_WM, combo2_down_WM, combo3_down_WM)
combo_down_WM <- unique(combo_down_WM[,1:3])

Decreasing <- rbind(combo_down_CA1, combo_down_CA3_4, combo_down_GCL, combo_down_ML, combo_down_SGZ,
    combo_down_SL, combo_down_SLM, combo_down_SR, combo_down_WM)

# Isolate rows repeated > 2 times for stringency in negative aging signature

Common_decreasing <- Decreasing %>%
  group_by_all() %>%
  filter(n() > 2) %>%
  as.data.frame() %>%
  unique()

# Create signed gene set for aging

# Check that there is no overlap
intersect(Common_increasing$gene_id,Common_decreasing$gene_id)

Common_increasing$sign <- 1
Common_decreasing$sign <- -1

CAS <- rbind(Common_increasing, Common_decreasing)

# Double check for no duplicates
CAS <- unique(CAS)

# directory to save aging gene set results
dir_outputs <- here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results")

fn_out4 <- file.path(dir_outputs, "CAS__BayesSpace_list")

# Export summary as .csv file
write.csv(CAS, fn_out4, row.names = FALSE)

# To pick up where you left off
# CAS <- read.csv(file = here("processed-data","pseudobulk_spe", "pseudoBulkDGE_results", "CAS__BayesSpace_list.csv"))

# Create a named vector or VISION gene signature creation

AgingSignature <- setNames(CAS$sign, CAS$gene_id)

########################
# VISION analysis of CAS
########################

sig <- createGeneSignature(name = "CommonAgingSignature", sigData = AgingSignature)

mySignatures <- c(sig)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# For unique labels use spe$key (some barcodes are duplicated as labels)
colnames(spe) <- spe$key

# Scale counts within a sample
n.umi <- colSums(assays(spe)$counts)

scaled_counts <- t(t(assays(spe)$counts) / n.umi) * median(n.umi)

vis <- Vision(scaled_counts,
              signatures = mySignatures)

# Set the number of threads when running parallel computations
# On Windows, this must either be omitted or set to 1
options(mc.cores = 1)

vis <- analyze(vis)

#The analysis has finished, and we'll remap now the scores for each of the signatures
  sigScores <- vis@SigScores[,1]
  sigRanks <- rank(sigScores)
  colData(spe)$CAS <- as.numeric(plyr::mapvalues(x = row.names(colData(spe)),
      from = names(sigScores), to = as.numeric(sigScores)))

# order spe observations according to age for better plotting
spe <- spe[, order(spe$age)]

# Plot CAS on tissue ordered by age
vis_grid_gene(
    spe = spe,
    geneid = "CAS",
    pdf = here::here("plots", "Trajectories", "CAS_BayesSpace_spotplot.pdf"),
    minCount = 0,
    viridis = FALSE,
    point_size = 1.5,
    spatial = TRUE,
    image_id = "lowres"
    )

# Setting up data for violin plots of CAS vs age_bin faceted by spatial domain

df <-
    data.frame(spe$key, spe$sample_id, spe$bayesSpace_harmony_10)

df <- df %>%
    mutate(
        BayesSpace = case_when(
            spe.bayesSpace_harmony_10 == 1 ~ "SLM",
            spe.bayesSpace_harmony_10 == 2 ~ "ML",
            spe.bayesSpace_harmony_10 == 3 ~ "CP",
            spe.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe.bayesSpace_harmony_10 == 5 ~ "SR",
            spe.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe.bayesSpace_harmony_10 == 7 ~ "GCL",
            spe.bayesSpace_harmony_10 == 8 ~ "SL",
            spe.bayesSpace_harmony_10 == 9 ~ "CA1",
            spe.bayesSpace_harmony_10 == 10 ~ "WM",
        )
    )

colData(spe)$BayesSpace <-
    factor(df$BayesSpace, levels = c("WM", "SLM", "SL", "ML", "CP", "SR", "SGZ",
        "CA3&4", "CA1", "GCL"))

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

saveRDS(spe,
    file = here::here("processed-data", "harmony_processed_spe", "harmony_CAS_BayesSpace_spe.rds"))

# Create datafame of cell proportions and DG layer
CAS_df <- as.data.frame(colData(spe)[, c(44:46)],
    row.names = spe$key)

# Add ages
CAS_df$age <- spe$age

# Remove CP cluster since that has higher variance and was not included in upstream DE analyses
# Remove Choroid Plexus cluster
CAS_df$spatialnum <- spe$bayesSpace_harmony_10
CAS_df = CAS_df[which(CAS_df$spatialnum != "3"), ]

strip <- strip_themed(background_x = elem_list_rect(fill = c("#2ED9FF", "#5A5156", "#90AD1C", "#E4E1E3",
    "#FE00FA", "#1CFFCE", "#FEAF16", "#16FF32", "#B00068")))

bay_colors <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SR" = "#FE00FA",
    "SGZ" = "#1CFFCE", "GCL" = "#B00068", "SL" = "#90AD1C", "CA1" =  "#16FF32", "WM" = "#2ED9FF")

pdf(file = here::here("plots", "Trajectories", "CAS_BayesSpace_vs_age.pdf"), width = 5, height = 5)

ggplot(CAS_df, aes(x = age, y = CAS)) +
    geom_boxplot(width = 10, notch = F, outlier.colour = NA,  mapping = aes(fill = as.factor(age))) +
    geom_smooth(data = CAS_df, mapping = aes(x = age, y = CAS), color='red3', se = F, method = 'loess') +
    labs(x = "age", y = "CAS (a.u.)") +
    facet_wrap2(~ BayesSpace, nrow = 3, ncol = 3, strip = strip) +
    theme(axis.text.x = element_text(size = 7), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank()) +
    guides(fill = FALSE) +
    ylim(-0.6, 2) +
    theme_classic()

ggplot(CAS_df, aes(x = CAS_df$age_bin, y = CAS_df$CAS)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width = 0.1) +
    geom_smooth(data = CAS_df, mapping = aes(x = as.numeric(age_bin), y = CAS), color='red3', se = F, method = 'loess') +
    labs(x = "age", y = "CAS (a.u.)") +
    ggtitle("Common Aging Signature") +
    theme_classic() +
    ylim(-0.6, 1) +
    facet_wrap2(~ BayesSpace, nrow = 3, ncol = 3, strip = strip) +
    theme(axis.text.x = element_text(size = 7), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(CAS_df, aes(x = age, y = CAS, color = BayesSpace)) +
    geom_smooth(data = CAS_df, mapping = aes(x = age, y = CAS), se = F,  method = 'loess') +
    labs(y = "CAS (a.u.)") +
    scale_color_manual(values = bay_colors) +
    scale_y_continuous(limits = c(-0.3, 0.8)) +
    theme_classic()

ggplot(CAS_df, aes(x = age, y = CAS, color = BayesSpace)) +
    geom_smooth(data = CAS_df, mapping = aes(x = age, y = CAS), se = F,  method = 'lm') +
    labs(y = "CAS (a.u.)") +
    scale_color_manual(values = bay_colors) +
    scale_y_continuous(limits = c(-0.3, 1.3)) +
    theme_classic()

dev.off()

# Linear modeling for CAS velocities

# set the spatial domain information to a new, character-type column, as the factors can confuse the linear models
CAS_df$spatial_domain <- as.character(CAS_df$BayesSpace)

# Calculate a linear model that allows for an interaction between age and spatial domain
m_interaction <- lm(CAS ~ age*spatial_domain, data = CAS_df)

# Setup a dataframe that contains the intercepts for each tissue for later inspection of the offset
coefficient_dataframe <- as.data.frame(m_interaction$coefficients)
colnames(coefficient_dataframe) <- 'coefficients'
coefficient_dataframe <- coefficient_dataframe[grep('^spatial_domain', row.names(coefficient_dataframe), value = T),,drop=F]
coefficient_dataframe$spatial_domain <- gsub('spatial_domain', '',row.names(coefficient_dataframe))
colnames(coefficient_dataframe)[1] <- 'offset'

# Use the function lstrends (from the lsmeans package) to allow domain-to-domain comparison
m_lst <- lstrends(m_interaction, "spatial_domain", var="age")

#Create pairwise comparisons between slopes and test using the pairs function and store that result
#This lets us know which slopes are statistically significantly different from one another
slope_pair_wiseTestResults <- as.data.frame(pairs(m_lst))

slope_pair_wiseTestResults

# contrast        estimate           SE    df t.ratio p.value
# CA1 - CA3&4  0.000281832 0.0001592537 66648   1.770  0.7024
# CA1 - GCL    0.001910985 0.0001751881 66648  10.908  <.0001
# CA1 - ML    -0.000102862 0.0001682034 66648  -0.612  0.9996
# CA1 - SGZ    0.000265943 0.0001715036 66648   1.551  0.8313
# CA1 - SL    -0.003131012 0.0001765065 66648 -17.739  <.0001
# CA1 - SLM   -0.002928358 0.0001597608 66648 -18.330  <.0001
# CA1 - SR    -0.000056485 0.0001701715 66648  -0.332  1.0000
# CA1 - WM    -0.002947631 0.0001644324 66648 -17.926  <.0001
# CA3&4 - GCL  0.001629153 0.0001496107 66648  10.889  <.0001
# CA3&4 - ML  -0.000384694 0.0001413678 66648  -2.721  0.1402
# CA3&4 - SGZ -0.000015889 0.0001452790 66648  -0.109  1.0000
# CA3&4 - SL  -0.003412843 0.0001511524 66648 -22.579  <.0001
# CA3&4 - SLM -0.003210190 0.0001312098 66648 -24.466  <.0001
# CA3&4 - SR  -0.000338317 0.0001437039 66648  -2.354  0.3096
# CA3&4 - WM  -0.003229463 0.0001368594 66648 -23.597  <.0001
# GCL - ML    -0.002013847 0.0001591037 66648 -12.657  <.0001
# GCL - SGZ   -0.001645041 0.0001625888 66648 -10.118  <.0001
# GCL - SL    -0.005041996 0.0001678576 66648 -30.037  <.0001
# GCL - SLM   -0.004839343 0.0001501504 66648 -32.230  <.0001
# GCL - SR    -0.001967470 0.0001611830 66648 -12.206  <.0001
# GCL - WM    -0.004858616 0.0001551116 66648 -31.323  <.0001
# ML - SGZ     0.000368806 0.0001550375 66648   2.379  0.2955
# ML - SL     -0.003028149 0.0001605543 66648 -18.861  <.0001
# ML - SLM    -0.002825496 0.0001419389 66648 -19.906  <.0001
# ML - SR      0.000046377 0.0001535625 66648   0.302  1.0000
# ML - WM     -0.002844769 0.0001471772 66648 -19.329  <.0001
# SGZ - SL    -0.003396955 0.0001640086 66648 -20.712  <.0001
# SGZ - SLM   -0.003194301 0.0001458348 66648 -21.904  <.0001
# SGZ - SR    -0.000322428 0.0001571706 66648  -2.051  0.5071
# SGZ - WM    -0.003213574 0.0001509379 66648 -21.291  <.0001
# SL - SLM     0.000202654 0.0001516866 66648   1.336  0.9207
# SL - SR      0.003074526 0.0001626150 66648  18.907  <.0001
# SL - WM      0.000183380 0.0001565992 66648   1.171  0.9626
# SLM - SR     0.002871873 0.0001442657 66648  19.907  <.0001
# SLM - WM    -0.000019273 0.0001374492 66648  -0.140  1.0000
# SR - WM     -0.002891146 0.0001494225 66648 -19.349  <.0001

#P value adjustment: tukey method for comparing a family of 9 estimates

# Store it in a separate dataframe to hold the information about the slopes
Visium_slope_dataframe <- as.data.frame(m_lst)

# Take the coefficient dataframe that we generate above together with the slope dataframe
Visium_coefficient_dataframe <- merge.data.frame(Visium_slope_dataframe, coefficient_dataframe, by = 'spatial_domain')

# Test correlation between Baseline CAS and slope
cor.test(Visium_coefficient_dataframe$offset, Visium_coefficient_dataframe$age.trend)

# Fit a linear model to the slope data for plot annotation
fit_trend <- lm(age.trend ~ offset, data = Visium_coefficient_dataframe)

# Extract the coefficients and R2 value
intercept_trend <- coef(fit_trend)[1]
slope_trend <- coef(fit_trend)[2]
r2_trend <- summary(fit_trend)$r.squared

pdf(file = here::here("plots", "Trajectories", "CAS_BayesSpace_slope_vs_baseline.pdf"), width = 7, height = 5)

ggplot(Visium_coefficient_dataframe, aes(x = offset, y = age.trend, color = spatial_domain)) +
    scale_color_manual(values = bay_colors) +
    geom_point(size = 10)  +
    labs(x = "Offset (a.u.)", y = "CAS Slope (a.u.)") +
    geom_smooth(method = 'lm', color='black', linetype = 'dashed', se = F) +
    scale_y_continuous(limits = c(0.002,0.01))  +
    scale_x_continuous(limits = c(-0.1, 0.7)) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black'),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 0.05, y = 0.0095,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_trend, r2_trend),
        size = 4) +
    annotate("text", x = 0.05, y = 0.0090,
        label = sprintf("R2 = %.2f", r2_trend), size = 4) +
    theme_classic()

dev.off()

# Compare CAS slope per spatial domain in barplot & add adjusted p-value (above) later

pdf(file = here::here("plots", "Trajectories", "CAS_BayesSpace_slope_spatial_domains.pdf"), width = 9, height = 5)

ggplot(Visium_slope_dataframe, aes(x = reorder(spatial_domain, -age.trend), y = age.trend, fill = age.trend)) +
    geom_bar(position = position_dodge(), stat = "identity", mapping = aes(fill = age.trend)) +
     xlab("") +
     ylab("CAS slope")  +
    labs(fill = "CAS slope") +
    scale_fill_gradient2(low = brewer.pal(11, "RdYlBu")[10], midpoint = quantile(Visium_coefficient_dataframe$age.trend)[3],
        mid = brewer.pal(11, "RdYlBu")[6], high = brewer.pal(11, "RdYlBu")[2])+
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = .2,                    # Width of the error bars
                  position  =position_dodge(.9)) +
     scale_y_continuous(expand = c(0,0)) +
     theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1)) +
    theme_classic()

dev.off()
