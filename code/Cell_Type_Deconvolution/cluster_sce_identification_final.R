######################################################
# spatial_DG_lifespan project
# Identifying clusters of sce object with marker genes
# Anthony Ramnauth, Nov 06 2022
######################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(dplyr)
    library(tidyr)
    library(forcats)
    library(ggplot2)
    library(RColorBrewer)
    library(ComplexHeatmap)
	library(sessioninfo)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_clustered.rds"))

colData(sce)$cell_type <- as.factor(colData(sce)$cell_type)

# number of nuclei per cluster and Dataset
table(colData(sce)$cell_type)
table(colData(sce)$cell_type, colData(sce)$Dataset)

########################################################
# Using find markers method from Tran et al. 2021 NEURON
# Modified since I don't have donor information.
########################################################

# Run pairwise t-tests
markers_t_pw <- findMarkers(sce, groups = sce$cell_type, block = sce$Dataset, test = "t",
    direction = "up")

sapply(markers_t_pw, function(x){table(x$FDR<0.05)})
#      Astro1 Astro2 Astro3 Astro4 Astro5 Astro6  COP1 Endo_Mur1 Endo_Mur2
#FALSE  49789  57704  55473  55103  58543  63870 64546     62030     58517
#TRUE   18969  11054  13285  13655  10215   4888  4212      6728     10241
#      ExctN1 ExctN2 ExctN3 ExctN4 ExctN5 ExctN6 ExctN7 ExctN8 ExctN9 Immun1
#FALSE  53555  47653  47188  47407  55763  51717  44520  48775  55595  59112
#TRUE   15203  21105  21570  21351  12995  17041  24238  19983  13163   9646
#      Immun2 Immun3 Immun4 Immun5 IntN1 IntN2 IntN3 IntN4 IntN5 IntN6 IntN7
#FALSE  58441  60606  65553  59494 59135 58964 58100 57907 58376 59808 62667
#TRUE   10317   8152   3205   9264  9623  9794 10658 10851 10382  8950  6091
#      Oligo1 Oligo2 Oligo3 Oligo4  OPC1
#FALSE  43559  52060  55995  60752 52388#
#TRUE   25199  16698  12763   8006 16370

# Run pairwise Wilcoxon rank sum tests
markers_wilcox_pw <- findMarkers(sce, groups = sce$cell_type, block = sce$Dataset, test = "wilcox",
    direction = "up")

sapply(markers_wilcox_pw, function(x){table(x$FDR<0.05)})
#      Astro1 Astro2 Astro3 Astro4 Astro5 Astro6  COP1 Endo_Mur1 Endo_Mur2
#FALSE  61056  66267  63384  59800  66721  64994 66172     64617     64237
#TRUE    7702   2491   5374   8958   2037   3764  2586      4141      4521
#      ExctN1 ExctN2 ExctN3 ExctN4 ExctN5 ExctN6 ExctN7 ExctN8 ExctN9 Immun1
#FALSE  58335  57868  56739  59890  58158  62736  57196  57573  64399  64032
#TRUE   10423  10890  12019   8868  10600   6022  11562  11185   4359   4726
#      Immun2 Immun3 Immun4 Immun5 IntN1 IntN2 IntN3 IntN4 IntN5 IntN6 IntN7
#FALSE  63492  63859  63190  64775 66253 66758 63340 64107 62949 66422 67030
#TRUE    5266   4899   5568   3983  2505  2000  5418  4651  5809  2336  1728
#      Oligo1 Oligo2 Oligo3 Oligo4  OPC1
#FALSE  60730  63122  66060  66719 61912
#TRUE    8028   5636   2698   2039  6846

# Wilcox test is more stringent with FDR values

# Run pairwise binomial tests
markers_binom_pw <- findMarkers(sce, groups = sce$cell_type, block = sce$Dataset, test = "binom",
    direction = "up")

sapply(markers_binom_pw, function(x){table(x$FDR<0.05)})
#      Astro1 Astro2 Astro3 Astro4 Astro5 Astro6  COP1 Endo_Mur1 Endo_Mur2
#FALSE  55376  61768  64018  57401  62609  65838 66374     61814     61529
#TRUE   13382   6990   4740  11357   6149   2920  2384      6944      7229
#      ExctN1 ExctN2 ExctN3 ExctN4 ExctN5 ExctN6 ExctN7 ExctN8 ExctN9 Immun1
#FALSE  56895  56810  55707  56291  56791  63182  52064  55125  55403  61621
#TRUE   11863  11948  13051  12467  11967   5576  16694  13633  13355   7137
#      Immun2 Immun3 Immun4 Immun5 IntN1 IntN2 IntN3 IntN4 IntN5 IntN6 IntN7
#FALSE  62503  63632  63966  63831 63664 66207 57880 62530 58009 66189 66935
#TRUE    6255   5126   4792   4927  5094  2551 10878  6228 10749  2569  1823
#      Oligo1 Oligo2 Oligo3 Oligo4  OPC1
#FALSE  56173  61368  62728  66365 58411
#TRUE   12585   7390   6030   2393 10347

#################################################
# Make .csv lists of top markers for each cluster
#################################################

# Make a data frame summary
Astro_1_wilcox_pw <- markers_wilcox_pw[[1]]
Astro_1_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Astro_1_wilcox_pw),
    rank = Astro_1_wilcox_pw$Top,
    p_val = Astro_1_wilcox_pw$p.value,
    FDR = Astro_1_wilcox_pw$FDR
)

Astro_1_wilcox_pw_summary <- Astro_1_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
dir_outputs <- here("processed-data", "sce", "cluster_markers")
fn_out_1 <- file.path(dir_outputs, "Astro_1_wilcox_pw_results")

# Export summary as .csv file
write.csv(Astro_1_wilcox_pw_summary,fn_out_1, row.names = FALSE)

# Make a data frame summary
Astro_2_wilcox_pw <- markers_wilcox_pw[[2]]
Astro_2_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Astro_2_wilcox_pw),
    rank = Astro_2_wilcox_pw$Top,
    p_val = Astro_2_wilcox_pw$p.value,
    FDR = Astro_2_wilcox_pw$FDR
)

Astro_2_wilcox_pw_summary <- Astro_2_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_2 <- file.path(dir_outputs, "Astro_2_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Astro_2_wilcox_pw_summary,fn_out_2, row.names = FALSE)

# Make a data frame summary
Astro_3_wilcox_pw <- markers_wilcox_pw[[3]]
Astro_3_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Astro_3_wilcox_pw),
    rank = Astro_3_wilcox_pw$Top,
    p_val = Astro_3_wilcox_pw$p.value,
    FDR = Astro_3_wilcox_pw$FDR
)

Astro_3_wilcox_pw_summary <- Astro_3_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_3 <- file.path(dir_outputs, "Astro_3_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Astro_3_wilcox_pw_summary,fn_out_3, row.names = FALSE)

# Make a data frame summary
Astro_4_wilcox_pw <- markers_wilcox_pw[[4]]
Astro_4_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Astro_4_wilcox_pw),
    rank = Astro_4_wilcox_pw$Top,
    p_val = Astro_4_wilcox_pw$p.value,
    FDR = Astro_4_wilcox_pw$FDR
)

Astro_4_wilcox_pw_summary <- Astro_4_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_4 <- file.path(dir_outputs, "Astro_4_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Astro_4_wilcox_pw_summary,fn_out_4, row.names = FALSE)

# Make a data frame summary
Astro_5_wilcox_pw <- markers_wilcox_pw[[5]]
Astro_5_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Astro_5_wilcox_pw),
    rank = Astro_5_wilcox_pw$Top,
    p_val = Astro_5_wilcox_pw$p.value,
    FDR = Astro_5_wilcox_pw$FDR
)

Astro_5_wilcox_pw_summary <- Astro_5_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_5 <- file.path(dir_outputs, "Astro_5_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Astro_5_wilcox_pw_summary,fn_out_5, row.names = FALSE)

# Make a data frame summary
Astro_6_wilcox_pw <- markers_wilcox_pw[[6]]
Astro_6_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Astro_6_wilcox_pw),
    rank = Astro_6_wilcox_pw$Top,
    p_val = Astro_6_wilcox_pw$p.value,
    FDR = Astro_6_wilcox_pw$FDR
)

Astro_6_wilcox_pw_summary <- Astro_6_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_6 <- file.path(dir_outputs, "Astro_6_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Astro_6_wilcox_pw_summary,fn_out_6, row.names = FALSE)

# Make a data frame summary
COP_wilcox_pw <- markers_wilcox_pw[[7]]
COP_wilcox_pw_summary <- data.frame(
    gene_name = rownames(COP_wilcox_pw),
    rank = COP_wilcox_pw$Top,
    p_val = COP_wilcox_pw$p.value,
    FDR = COP_wilcox_pw$FDR
)

COP_wilcox_pw_summary <- COP_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_7 <- file.path(dir_outputs, "COP_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(COP_wilcox_pw_summary,fn_out_7, row.names = FALSE)

# Make a data frame summary
Endo_Mur1_wilcox_pw <- markers_wilcox_pw[[8]]
Endo_Mur1_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Endo_Mur1_wilcox_pw),
    rank = Endo_Mur1_wilcox_pw$Top,
    p_val = Endo_Mur1_wilcox_pw$p.value,
    FDR = Endo_Mur1_wilcox_pw$FDR
)

Endo_Mur1_wilcox_summary <- Endo_Mur1_wilcox_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_8 <- file.path(dir_outputs, "Endo_Mur1_wilcox_test_results")

# Export summary as .csv file
write.csv(Endo_Mur1_wilcox_pw_summary,fn_out_8, row.names = FALSE)

# Make a data frame summary
Endo_Mur2_wilcox_pw <- markers_wilcox_pw[[9]]
Endo_Mur2_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Endo_Mur2_wilcox_pw),
    rank = Endo_Mur2_wilcox_pw$Top,
    p_val = Endo_Mur2_wilcox_pw$p.value,
    FDR = Endo_Mur2_wilcox_pw$FDR
)

Endo_Mur2_wilcox_pw_summary <- Endo_Mur2_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_9 <- file.path(dir_outputs, "Endo_Mur2_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Endo_Mur2_wilcox_pw_summary,fn_out_9, row.names = FALSE)

# Make a data frame summary
ExctN1_wilcox_pw <- markers_wilcox_pw[[10]]
ExctN1_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN1_wilcox_pw),
    rank = ExctN1_wilcox_pw$Top,
    p_val = ExctN1_wilcox_pw$p.value,
    FDR = ExctN1_wilcox_pw$FDR
)

ExctN1_wilcox_pw_summary <- ExctN1_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_10 <- file.path(dir_outputs, "ExctN1_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN1_wilcox_pw_summary,fn_out_10, row.names = FALSE)

# Make a data frame summary
ExctN2_wilcox_pw <- markers_wilcox_pw[[11]]
ExctN2_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN2_wilcox_pw),
    rank = ExctN2_wilcox_pw$Top,
    p_val = ExctN2_wilcox_pw$p.value,
    FDR = ExctN2_wilcox_pw$FDR
)

ExctN2_wilcox_pw_summary <- ExctN2_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_11 <- file.path(dir_outputs, "ExctN2_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN2_wilcox_pw_summary,fn_out_11, row.names = FALSE)

# Make a data frame summary
ExctN3_wilcox_pw <- markers_wilcox_pw[[12]]
ExctN3_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN3_wilcox_pw),
    rank = ExctN3_wilcox_pw$Top,
    p_val = ExctN3_wilcox_pw$p.value,
    FDR = ExctN3_wilcox_pw$FDR
)

ExctN3_wilcox_pw_summary <- ExctN3_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_12 <- file.path(dir_outputs, "ExctN3_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN3_wilcox_pw_summary,fn_out_12, row.names = FALSE)

# Make a data frame summary
ExctN4_wilcox_pw <- markers_wilcox_pw[[13]]
ExctN4_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN4_wilcox_pw),
    rank = ExctN4_wilcox_pw$Top,
    p_val = ExctN4_wilcox_pw$p.value,
    FDR = ExctN4_wilcox_pw$FDR
)

ExctN4_wilcox_pw_summary <- ExctN4_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_13 <- file.path(dir_outputs, "ExctN4_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN4_wilcox_pw_summary,fn_out_13, row.names = FALSE)

# Make a data frame summary
ExctN5_wilcox_pw <- markers_wilcox_pw[[14]]
ExctN5_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN5_wilcox_pw),
    rank = ExctN5_wilcox_pw$Top,
    p_val = ExctN5_wilcox_pw$p.value,
    FDR = ExctN5_wilcox_pw$FDR
)

ExctN5_wilcox_pw_summary <- ExctN5_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_14 <- file.path(dir_outputs, "ExctN5_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN5_wilcox_pw_summary,fn_out_14, row.names = FALSE)

# Make a data frame summary
ExctN6_wilcox_pw <- markers_wilcox_pw[[15]]
ExctN6_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN6_wilcox_pw),
    rank = ExctN6_wilcox_pw$Top,
    p_val = ExctN6_wilcox_pw$p.value,
    FDR = ExctN6_wilcox_pw$FDR
)

ExctN6_wilcox_pw_summary <- ExctN6_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_15 <- file.path(dir_outputs, "ExctN6_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN6_wilcox_pw_summary,fn_out_15, row.names = FALSE)

# Make a data frame summary
ExctN7_wilcox_pw <- markers_wilcox_pw[[16]]
ExctN7_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN7_wilcox_pw),
    rank = ExctN7_wilcox_pw$Top,
    p_val = ExctN7_wilcox_pw$p.value,
    FDR = ExctN7_wilcox_pw$FDR
)

ExctN7_wilcox_pw_summary <- ExctN7_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_16 <- file.path(dir_outputs, "ExctN7_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN7_wilcox_pw_summary,fn_out_16, row.names = FALSE)

# Make a data frame summary
ExctN8_wilcox_pw <- markers_wilcox_pw[[17]]
ExctN8_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN8_wilcox_pw),
    rank = ExctN8_wilcox_pw$Top,
    p_val = ExctN8_wilcox_pw$p.value,
    FDR = ExctN8_wilcox_pw$FDR
)

ExctN8_wilcox_pw_summary <- ExctN8_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_17 <- file.path(dir_outputs, "ExctN8_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN8_wilcox_pw_summary,fn_out_17, row.names = FALSE)

# Make a data frame summary
ExctN9_wilcox_pw <- markers_wilcox_pw[[18]]
ExctN9_wilcox_pw_summary <- data.frame(
    gene_name = rownames(ExctN9_wilcox_pw),
    rank = ExctN9_wilcox_pw$Top,
    p_val = ExctN9_wilcox_pw$p.value,
    FDR = ExctN9_wilcox_pw$FDR
)

ExctN9_wilcox_pw_summary <- ExctN9_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_18 <- file.path(dir_outputs, "ExctN9_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(ExctN9_wilcox_pw_summary,fn_out_18, row.names = FALSE)

# Make a data frame summary
Immun1_wilcox_pw <- markers_wilcox_pw[[19]]
Immun1_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Immun1_wilcox_pw),
    rank = Immun1_wilcox_pw$Top,
    p_val = Immun1_wilcox_pw$p.value,
    FDR = Immun1_wilcox_pw$FDR
)

Immun1_wilcox_pw_summary <- Immun1_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_19 <- file.path(dir_outputs, "Immun1_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Immun1_wilcox_pw_summary,fn_out_19, row.names = FALSE)

# Make a data frame summary
Immun2_wilcox_pw <- markers_wilcox_pw[[20]]
Immun2_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Immun2_wilcox_pw),
    rank = Immun2_wilcox_pw$Top,
    p_val = Immun2_wilcox_pw$p.value,
    FDR = Immun2_wilcox_pw$FDR
)

Immun2_wilcox_pw_summary <- Immun2_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_20 <- file.path(dir_outputs, "Immun2_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Immun2_wilcox_pw_summary,fn_out_20, row.names = FALSE)

# Make a data frame summary
Immun3_wilcox_pw <- markers_wilcox_pw[[21]]
Immun3_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Immun3_wilcox_pw),
    rank = Immun3_wilcox_pw$Top,
    p_val = Immun3_wilcox_pw$p.value,
    FDR = Immun3_wilcox_pw$FDR
)

Immun3_wilcox_pw_summary <- Immun3_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_21 <- file.path(dir_outputs, "Immun3_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Immun3_wilcox_pw_summary,fn_out_21, row.names = FALSE)

# Make a data frame summary
Immun4_wilcox_pw <- markers_wilcox_pw[[22]]
Immun4_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Immun4_wilcox_pw),
    rank = Immun4_wilcox_pw$Top,
    p_val = Immun4_wilcox_pw$p.value,
    FDR = Immun4_wilcox_pw$FDR
)

Immun4_wilcox_pw_summary <- Immun4_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_22 <- file.path(dir_outputs, "Immun4_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Immun4_wilcox_pw_summary,fn_out_22, row.names = FALSE)

# Make a data frame summary
Immun5_wilcox_pw <- markers_wilcox_pw[[23]]
Immun5_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Immun5_wilcox_pw),
    rank = Immun5_wilcox_pw$Top,
    p_val = Immun5_wilcox_pw$p.value,
    FDR = Immun5_wilcox_pw$FDR
)

Immun5_wilcox_pw_summary <- Immun5_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_23 <- file.path(dir_outputs, "Immun5_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Immun5_wilcox_pw_summary,fn_out_23, row.names = FALSE)

# Make a data frame summary
IntN1_wilcox_pw <- markers_wilcox_pw[[24]]
IntN1_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN1_wilcox_pw),
    rank = IntN1_wilcox_pw$Top,
    p_val = IntN1_wilcox_pw$p.value,
    FDR = IntN1_wilcox_pw$FDR
)

IntN1_wilcox_pw_summary <- IntN1_wilcox_pw_summary %>%
    filter(p_val < 0.2)

# directory to save lists
fn_out_24 <- file.path(dir_outputs, "IntN1_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN1_wilcox_pw_summary,fn_out_24, row.names = FALSE)

# Make a data frame summary
IntN2_wilcox_pw <- markers_wilcox_pw[[25]]
IntN2_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN2_wilcox_pw),
    rank = IntN2_wilcox_pw$Top,
    p_val = IntN2_wilcox_pw$p.value,
    FDR = IntN2_wilcox_pw$FDR
)

IntN2_wilcox_pw_summary <- IntN2_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_25 <- file.path(dir_outputs, "IntN2_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN2_wilcox_pw_summary,fn_out_25, row.names = FALSE)

# Make a data frame summary
IntN3_wilcox_pw <- markers_wilcox_pw[[26]]
IntN3_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN3_wilcox_pw),
    rank = IntN3_wilcox_pw$Top,
    p_val = IntN3_wilcox_pw$p.value,
    FDR = IntN3_wilcox_pw$FDR
)

IntN3_wilcox_pw_summary <- IntN3_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_26 <- file.path(dir_outputs, "IntN3_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN3_wilcox_pw_summary,fn_out_26, row.names = FALSE)

# Make a data frame summary
IntN4_wilcox_pw <- markers_wilcox_pw[[27]]
IntN4_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN4_wilcox_pw),
    rank = IntN4_wilcox_pw$Top,
    p_val = IntN4_wilcox_pw$p.value,
    FDR = IntN4_wilcox_pw$FDR
)

IntN4_wilcox_pw_summary <- IntN4_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_27 <- file.path(dir_outputs, "IntN4_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN4_wilcox_pw_summary,fn_out_27, row.names = FALSE)

# Make a data frame summary
IntN5_wilcox_pw <- markers_wilcox_pw[[28]]
IntN5_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN5_wilcox_pw),
    rank = IntN5_wilcox_pw$Top,
    p_val = IntN5_wilcox_pw$p.value,
    FDR = IntN5_wilcox_pw$FDR
)

IntN5_wilcox_pw_summary <- IntN5_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_28 <- file.path(dir_outputs, "IntN5_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN5_wilcox_pw_summary,fn_out_28, row.names = FALSE)

# Make a data frame summary
IntN6_wilcox_pw <- markers_wilcox_pw[[29]]
IntN6_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN6_wilcox_pw),
    rank = IntN6_wilcox_pw$Top,
    p_val = IntN6_wilcox_pw$p.value,
    FDR = IntN6_wilcox_pw$FDR
)

IntN6_wilcox_pw_summary <- IntN6_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_29 <- file.path(dir_outputs, "IntN6_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN6_wilcox_pw_summary,fn_out_29, row.names = FALSE)

# Make a data frame summary
IntN7_wilcox_pw <- markers_wilcox_pw[[30]]
IntN7_wilcox_pw_summary <- data.frame(
    gene_name = rownames(IntN7_wilcox_pw),
    rank = IntN7_wilcox_pw$Top,
    p_val = IntN7_wilcox_pw$p.value,
    FDR = IntN7_wilcox_pw$FDR
)

IntN7_wilcox_pw_summary <- IntN7_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_30 <- file.path(dir_outputs, "IntN7_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(IntN7_wilcox_pw_summary,fn_out_30, row.names = FALSE)

# Make a data frame summary
Oligo1_wilcox_pw <- markers_wilcox_pw[[31]]
Oligo1_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Oligo1_wilcox_pw),
    rank = Oligo1_wilcox_pw$Top,
    p_val = Oligo1_wilcox_pw$p.value,
    FDR = Oligo1_wilcox_pw$FDR
)

Oligo1_wilcox_pw_summary <- Oligo1_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_31 <- file.path(dir_outputs, "Oligo1_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Oligo1_wilcox_pw_summary,fn_out_31, row.names = FALSE)

# Make a data frame summary
Oligo2_wilcox_pw <- markers_wilcox_pw[[32]]
Oligo2_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Oligo2_wilcox_pw),
    rank = Oligo2_wilcox_pw$Top,
    p_val = Oligo2_wilcox_pw$p.value,
    FDR = Oligo2_wilcox_pw$FDR
)

Oligo2_wilcox_pw_summary <- Oligo2_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_32 <- file.path(dir_outputs, "Oligo2_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Oligo2_wilcox_pw_summary,fn_out_32, row.names = FALSE)

# Make a data frame summary
Oligo3_wilcox_pw <- markers_wilcox_pw[[33]]
Oligo3_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Oligo3_wilcox_pw),
    rank = Oligo3_wilcox_pw$Top,
    p_val = Oligo3_wilcox_pw$p.value,
    FDR = Oligo3_wilcox_pw$FDR
)

Oligo3_wilcox_pw_summary <- Oligo3_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_33 <- file.path(dir_outputs, "Oligo3_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Oligo3_wilcox_pw_summary,fn_out_33, row.names = FALSE)

# Make a data frame summary
Oligo4_wilcox_pw <- markers_wilcox_pw[[34]]
Oligo4_wilcox_pw_summary <- data.frame(
    gene_name = rownames(Oligo4_wilcox_pw),
    rank = Oligo4_wilcox_pw$Top,
    p_val = Oligo4_wilcox_pw$p.value,
    FDR = Oligo4_wilcox_pw$FDR
)

Oligo4_wilcox_pw_summary <- Oligo4_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_34 <- file.path(dir_outputs, "Oligo4_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(Oligo4_wilcox_pw_summary,fn_out_34, row.names = FALSE)

# Make a data frame summary
OPC_wilcox_pw <- markers_wilcox_pw[[35]]
OPC_wilcox_pw_summary <- data.frame(
    gene_name = rownames(OPC_wilcox_pw),
    rank = OPC_wilcox_pw$Top,
    p_val = OPC_wilcox_pw$p.value,
    FDR = OPC_wilcox_pw$FDR
)

OPC_wilcox_pw_summary <- OPC_wilcox_pw_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

# directory to save lists
fn_out_35 <- file.path(dir_outputs, "OPC_wilcox_pw_test_results")

# Export summary as .csv file
write.csv(OPC_wilcox_pw_summary,fn_out_35, row.names = FALSE)

# ------------
# Marker genes
# ------------

# Markers chosen from the publications of each dataset

markers <- c(
    ## Neurons
    "RBFOX3", "SNAP25", "SYT1",
    ## Excitatory Neurons
    "SLC17A7",
    # GC
    "PROX1",
    # im GC
    "DCX", "BHLHE22", "STMN1",
    # Mossy Cells
    "ARHGAP24", "DLC1",
    # CA3 PNs
    "CFAP299", "SYN3",
    # CA2 PNs
    "HGF",
    # CA1 PNs
    "ACVR1C", "SYT13",
    # Sub PNs
    "ROBO1", "COL5A2",
    # Cajal?Retzius
    "RELN",
    ## Inhibitory Neurons
    "GAD1", "GAD2",
    # inhibitory subpopulations (Some from Lukas LC & Keri Martinowich 2022-07-22)
    "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", "PVALB", "CORT", "VIP", "NPY",
    "CRHBP", "CCK", "HTR3A", "NR2F2", "LAMP5",
    # Astrocytes
    "AQP4", "GFAP", "CHRDL1",
    # Oligodendrocytes
    "MOBP",
    # macrophages / microglia
    "CD163", "C3", "PTPRC", "C1QB",
    # OPCs
    "PDGFRA", "VCAN",
    # COP
    "GPR17", "ADAM33",
    # endothelial / mural (RBPMS)
    "CLDN5", "FLT1", "RBPMS",
    # T cells
    "SKAP1", "CD247",
    # Progenitors
    "PAX6", "HOPX", "EOMES"
)

# marker labels
marker_labels <- c(
  rep("Neuron", 3),
  rep("Excitatory", 1),
  rep("Granular_cells", 1),
  rep("IM_Granular_cells", 3),
  rep("Mossy_Cells", 2),
  rep("CA3", 2),
  rep("CA2", 1),
  rep("CA1", 2),
  rep("Sub", 2),
  rep("Cajal?Retzius", 1),
  rep("Inhibitory", 17),
  rep("Astrocytes", 3),
  rep("Oligodendrocytes", 1),
  rep("Macrophages_Microglia", 4),
  rep("OPCs", 2),
  rep("COP", 2),
  rep("Endothelial_Mural", 3),
  rep("T_cells", 2),
  rep("Progenitors", 3)
)

# colors: selected from tableau20 and tableau10medium
colors_markers <- list(marker = c(
  Neuron = "blue",
  Excitatory = "darkblue",
  Granular_cells = "blue1",
  IM_Granular_cells = "blue2",
  Mossy_Cells = "blue4",
  CA3 = "dodgerblue",
  CA2 = "dodgerblue3",
  CA1 = "blue3",
  Sub = "blue4",
  'Cajal?Retzius' = "gray47",
  Inhibitory = "green",
  Astrocytes = "yellow",
  Oligodendrocytes = "plum3",
  Macrophages_Microglia = "tan",
  OPCs = "goldenrod",
  COP = "goldenrod4",
  Endothelial_Mural = "red3",
  T_cells = "tan3",
  Progenitors = "cyan"
    ))

marker_labels <-
  factor(marker_labels, levels = unique(marker_labels))

# number of nuclei per cluster
n <- table(colLabels(sce))

# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$label)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$SYMBOL

hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers, i]))))

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

hm_mat <- t(scale_rows(t(hm_mat)))

# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels,
  show_annotation_name = FALSE,
  show_legend = TRUE,
  col = colors_markers
  )

pdf(file = here::here("plots", "sce_plots", "Heatmap_markers_sce.pdf"), width = 12, height = 8)

Heatmap(
  hm_mat,
  name = "z-score",
  column_title = "DG clusters mean marker expression",
  column_title_gp = gpar(fontface = "bold"),
  bottom_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_split = marker_labels,
  column_names_gp = gpar(fontface = "italic"),
  rect_gp = gpar(col = "gray50", lwd = 0.5)
    )

dev.off()

######
# UMAP
######

# cluster labels
cluster_pops <- list(
  ExctN = c(4, 9, 1, 6, 8, 10),
  InhbN = c(19, 7, 16, 11, 20),
  Astro = c(13),
  Oligo = c(15, 12, 2, 5, 18, 14),
  Macro_Micro_T = c(17),
  OPCs = c(3),
  Endo_Mural = c(21)
    )

label_merged <- fct_collapse(colData(sce)$label,
  ExctN = as.character(cluster_pops[[1]]),
  InhbN = as.character(cluster_pops[[2]]),
  Astro = as.character(cluster_pops[[3]]),
  Oligo = as.character(cluster_pops[[4]]),
  Macro_Micro_T = as.character(cluster_pops[[5]]),
  OPCs = as.character(cluster_pops[[6]]),
  Endo_Mural = as.character(cluster_pops[[7]])
    )

label_merged <- fct_relevel(label_merged,
  c("ExctN", "InhbN", "Astro", "Oligo", "Macro_Micro_T",
      "OPCs", "Endo_Mural"))

table(label_merged)
#label_merged
#        ExctN         InhbN         Astro         Oligo Macro_Micro_T
#        60373         24090         27609         69768         13613
#         OPCs    Endo_Mural
#        12732          1090

colData(sce)$label_merged <- label_merged

colors_clusters <- list(population = c(
  ExctN = "blue",
  InhbN = "green",
  Astro = "yellow",
  Oligo = "plum3",
  Macro_Micro_T = "tan",
  OPCs = "goldenrod",
  Endo_Mural = "red3")
    )

pdf(file = here::here("plots", "sce_plots", "Merged_cluster_plot_sce.pdf"))

plotReducedDim(sce, dimred = "UMAP.HARMONY", colour_by = "label_merged") +
  scale_color_manual(values = colors_clusters[[1]], name = "clusters (merged)") +
  theme_classic() +
  ggtitle("Combined HPC snRNAseq datasets clustering")

dev.off()

saveRDS(sce, file = here::here("processed-data", "sce", "sce_clustered.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
