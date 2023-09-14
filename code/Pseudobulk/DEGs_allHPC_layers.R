##################################################
# spatial_DG_lifespan project
# DEGs for all HPC BayesSpace clusters and age_bin
# Anthony Ramnauth, Aug 04 2023
##################################################

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(ggplot2)
})

# Make Data frame for age group, cluster, total, DEG
# # of significant DEGs taken from .csv files output from pseudoBulkDGE.R script

spatial_domain <- c("SLM", "ML", "CA3&4", "SR", "SGZ", "GCL", "SL", "CA1", "WM")


df <- data.frame(
    spatial_domain = rep(spatial_domain, 4),
    age = c(rep("infant", 9), rep("teen", 9), rep("adult", 9), rep("elderly", 9)),
    DEGs = c(863, 860, 280, 429, 302, 1578, 576, 806, 1318,
        0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 2, 0, 0, 4, 19, 1, 6, 0,
        53, 23, 15, 7, 15, 12, 15, 1, 8)
)

df$spatial_domain <- factor(df$spatial_domain, levels = c("SLM", "ML", "CA3&4", "SR", "SGZ", "GCL", "SL", "CA1", "WM"))
df$age <- factor(df$age, levels = c("infant", "teen", "adult", "elderly"))

pdf(file = here::here("plots", "pseudobulked", "Barplot_DEGs_allHPC_layers.pdf"))

ggplot(df, aes(x = spatial_domain, y = DEGs, fill = age)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("# Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

ggplot(df, aes(x = spatial_domain, y = DEGs, fill = age)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("# Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

dev.off()

