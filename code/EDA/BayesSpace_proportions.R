###############################
# spatial_DG_lifespan project
# Plotting cluster proportions
# Anthony Ramnauth, Aug 25 2023
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(spatialLIBD)
    library(dplyr)
})

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# order spe observations according to age
spe <- spe[, order(spe$age)]

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

df <- as.data.frame(colData(spe))

proportions_data_per_sample <- df %>%
  group_by(sample_id) %>%
  arrange(age) %>%
  count(bayesSpace_harmony_10) %>%
  mutate(prop = prop.table(n))

proportions_data_per_age <- df %>%
  group_by(age_bin) %>%
  count(bayesSpace_harmony_10) %>%
  mutate(prop = prop.table(n))

sample_prop <- as.data.frame(proportions_data_per_sample)
sample_prop$age <- c(rep(15.16, 10), rep(17.94, 10), rep(48.20, 10), rep(47.50, 10),
    rep(73.90, 10), rep(92.25, 10), rep(76.38, 10), rep(1.84, 10), rep(18.41, 10),
    rep(33.40, 10), rep(17.42, 10), rep(0.31, 10), rep(0.58, 10), rep(37.30, 10),
    rep(1.05, 10), rep(0.18, 10))

sample_prop <- sample_prop %>%
  arrange(age)

level_order <- c("Br8700", "Br8195", "Br8533", "Br8686", "Br6129_new", "Br1412", "Br8181",
    "Br2706", "Br6299_new", "Br6522", "Br8667", "Br3942", "Br2720", "Br5242", "Br6023", "Br5699_new")

age_prop <- as.data.frame(proportions_data_per_age)

pdf(file = here::here("plots", "BayesSpace_plots", "BayesSpace_spot_proportions.pdf"))

ggplot(sample_prop, aes(x = factor(sample_id, level = level_order), y = prop, fill = bayesSpace_harmony_10)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual("BayesSpace", values = c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")) +
    ylab("Proportion of Visium Spots") +
    xlab("sample ID") +
    theme(text = element_text(size = 14), axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))

ggplot(age_prop, aes(x = age_bin, y = prop, fill = bayesSpace_harmony_10)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual("BayesSpace", values = c("1" = "#5A5156", "2" = "#E4E1E3", "3" = "#DEA0FD", "4" = "#FEAF16",
    "5" = "#FE00FA", "6" = "#1CFFCE", "7" = "#B00068", "8" = "#90AD1C", "9" = "#16FF32",
    "10" = "#2ED9FF")) +
    ylab("Proportion of Visium Spots") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

dev.off()
