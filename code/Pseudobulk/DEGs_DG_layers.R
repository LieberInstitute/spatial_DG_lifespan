##############################################
# spatial_DG_lifespan project
# DEGs for each DG layer and age_bin
# Anthony Ramnauth, Aug 04 2023
##############################################

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(ggplot2)
    library(ggVennDiagram)
})

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

pdf(file = here::here("plots", "pseudobulked", "BLAH_Barplot_DEGs_DG_layers.pdf"))

ggplot(df, aes(x = layer, y = DEGs, fill = age)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_y_log10() +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("log10 number of Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

ggplot(df, aes(x = layer, y = DEGs, fill = age)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_y_log10() +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("log10 number of Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

dev.off()
