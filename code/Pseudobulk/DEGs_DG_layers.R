##############################################
# spatial_DG_lifespan project
# DEGs for each DG layer and age_bin
# Anthony Ramnauth, Aug 04 2023
##############################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(sessioninfo)
})

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Make Data frame for age group, cluster, total, DEG
# # of significant DEGs taken from .csv files output from pseudoBulkDGE.R script

layer <- c("ML", "GCL", "SGZ", "CA3&4")


df <- data.frame(
    layer = rep(layer, 4),
    age = c(rep("infant", 4), rep("teen", 4), rep("adult", 4), rep("elderly", 4)),
    DEGs = c(860, 1578, 302, 280, 0, 0, 1, 0, 2, 19, 4, 0, 25, 12, 15, 15)
)

df$layer <- factor(df$layer, levels = c("ML", "GCL", "SGZ", "CA3&4"))
df$age <- factor(df$age, levels = c("infant", "teen", "adult", "elderly"))

df

#   layer     age DEGs
#1     ML  infant  860
#2    GCL  infant 1578
#3    SGZ  infant  302
#4  CA3&4  infant  280
#5     ML    teen    0
#6    GCL    teen    0
#7    SGZ    teen    1
#8  CA3&4    teen    0
#9     ML   adult    2
#10   GCL   adult   19
#11   SGZ   adult    4
#12 CA3&4   adult    0
#13    ML elderly   25
#14   GCL elderly   12
#15   SGZ elderly   15
#16 CA3&4 elderly   15

pdf(file = here::here("plots", "pseudobulked", "Barplot_DEGs_DG_layers.pdf"), width = 5, height = 5)

ggplot(df, aes(x = layer, y = DEGs, fill = age)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_y_log10() +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("log10 number of Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    theme_classic()

ggplot(df, aes(x = layer, y = DEGs, fill = age)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_y_log10() +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("log10 number of Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    theme_classic()

dev.off()

########################################################################################################################

################################################################
# DEGs vs. # of Visium spots for each spatial domain + age group
################################################################

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

table(spe$bayesSpace_harmony_10, spe$age_bin)

#     Infant Teen Adult Elderly
#  1    2999 2251  2602    1716
#  2    2105 2019  1975    1176
#  3    1141  327   402     149
#  4    2528 3187  3109    1828
#  5    2964 2057  1784    1118
#  6    1686 1566  2009    1184
#  7    1260 1348  1737    1117
#  8    2827 1541  1595     892
#  9    1471 1267  1274    1109
#  10   2431 1816  1583    1535

df$spots <- c(2105, 1260, 1686, 2528, 2019, 1348, 1566, 3187, 1975, 1737, 2009, 3109, 1176, 1117, 1184, 1828)

# Fit a linear model to the log-transformed data
fit_spots <- lm(DEGs ~ spots, data = df)

# Extract the coefficients and R2 value
intercept_spots <- coef(fit_spots)[1]
slope_spots <- coef(fit_spots)[2]
r2_spots <- summary(fit_spots)$r.squared

###########################################################################################

#######################################################
# DEGs vs. # of UMI for each spatial domain + age group
#######################################################

coldata_spe <- as.data.frame(colData(spe))

# Convert column C to numeric
coldata_spe$new_age_bin <- as.numeric(coldata_spe$age_bin)

# Calculate the sum of UMIs limited by BayesSpace cluster and age_bin
sum_UMI <- c(
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 4]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 4]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 4]),
    sum(coldata_spe$sum_umi[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 4])
    )

df$UMI <- sum_UMI

# Fit a linear model to the log-transformed data
fit_UMI <- lm(DEGs ~ UMI, data = df)

# Extract the coefficients and R2 value
intercept_UMI <- coef(fit_UMI)[1]
slope_UMI <- coef(fit_UMI)[2]
r2_UMI <- summary(fit_UMI)$r.squared

#############################################################################################################

##################################################################
# DEGs vs. # of genes detected for each spatial domain + age group
##################################################################

sum_gene <- c(
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 1]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 2]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 3]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 2 & coldata_spe$new_age_bin == 4]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 7 & coldata_spe$new_age_bin == 4]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 6 & coldata_spe$new_age_bin == 4]),
    sum(coldata_spe$sum_gene[coldata_spe$bayesSpace_harmony_10 == 4 & coldata_spe$new_age_bin == 4])
    )

df$genes <- sum_gene

# Fit a linear model to the log-transformed data
fit_gene <- lm(DEGs ~ genes, data = df)

# Extract the coefficients and R2 value
intercept_gene <- coef(fit_gene)[1]
slope_gene <- coef(fit_gene)[2]
r2_gene <- summary(fit_gene)$r.squared

############################################################################################################

# Plot these as scatter plots

custom_colors <- c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")

pdf(file = here::here("plots", "pseudobulked", "DEGsvslayerandage.pdf"), width = 7, height = 5)

ggplot(df, aes(x = spots, y = DEGs, color = age)) +
    geom_point(size = 5) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_manual(values = custom_colors) +
    geom_text_repel(aes(label = layer),max.overlaps = 15) +
    labs(title = "DEGs vs spots across spatial domain + age group",
        x = "Number of spots (log)",
        y = "Number of significant DEGs (log)") +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 2500, y = 3000,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_spots, r2_spots),
        size = 4) +
    annotate("text", x = 2500, y = 2000,
        label = sprintf("R2 = %.2f", r2_spots), size = 4)

ggplot(df, aes(x = UMI, y = DEGs, color = age)) +
    geom_point(size = 5) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_manual(values = custom_colors) +
    geom_text_repel(aes(label = layer),max.overlaps = 15) +
    labs(title = "DEGs vs UMIs across spatial domain + age group",
        x = "Number of UMIs (log)",
        y = "Number of significant DEGs (log)") +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 2e+07, y = 3000,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_UMI, r2_UMI),
        size = 4) +
    annotate("text", x = 2e+07, y = 2000,
        label = sprintf("R2 = %.2f", r2_UMI), size = 4)

ggplot(df, aes(x = genes, y = DEGs, color = age)) +
    geom_point(size = 5) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_manual(values = custom_colors) +
    geom_text_repel(aes(label = layer),max.overlaps = 15) +
    labs(title = "DEGs vs genes detected across spatial domain + age group",
        x = "Number of genes (log)",
        y = "Number of significant DEGs (log)") +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size =10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 7.5e+06, y = 3000,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_gene, r2_gene),
        size = 4) +
    annotate("text", x = 7.5e+06, y = 2000,
        label = sprintf("R2 = %.2f", r2_gene), size = 4)

dev.off()
