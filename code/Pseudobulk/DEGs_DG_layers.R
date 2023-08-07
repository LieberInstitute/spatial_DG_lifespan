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
    library(sessioninfo)
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

pdf(file = here::here("plots", "pseudobulked", "Barplot_DEGs_DG_layers.pdf"))

ggplot(df, aes(x = layer, y = DEGs, fill = age)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("# Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

ggplot(df, aes(x = layer, y = DEGs, fill = age)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual("legend", values = c("infant" = "purple", "teen" = "blue", "adult" = "red", "elderly" = "forestgreen")) +
    ylab("# Differentially Expressed Genes") +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_classic()

dev.off()

########################################################################################################################

# Check overlap of ML and GCL DEGs (Are most just mRNA trafficking from GCL?)

# Load .csv files with the list of DE genes

Infant_ML_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "InfantvsNonInfant_BayesSpace2_DE.csv"))

Infant_GCL_DE_age_results <- read.csv(file = here::here("processed-data", "pseudobulk_spe",
    "pseudoBulkDGE_results", "InfantvsNonInfant_BayesSpace7_DE.csv"))

## List for Venn Diagram
infant_ML_GCL <- list(infant_ML = Infant_ML_DE_age_results$gene_name,
                   infant_GCL = Infant_GCL_DE_age_results$gene_name
    )

## Plot Venn Diagram

pdf(file = here::here("plots", "pseudobulked", "Venn_infant_ML_GCL.pdf"))

ggVennDiagram(infant_ML_GCL, label_alpha = 0, set_size = 7, label_size = 6) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  ggplot2::scale_fill_gradient(low = "white",high = "red")

dev.off()

# Get gene list non-overlapping DEGs for infant ML

infant_ML_only <- setdiff(Infant_ML_DE_age_results$gene_name, Infant_GCL_DE_age_results$gene_name)

# directory to save non-overlapping DEGs for infant ML
dir_outputs <- here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results")
fn_out <- file.path(dir_outputs, "infant_ML_notGCL_DEGs")

# Export summary as .csv file
write.csv(infant_ML_only, fn_out, row.names = FALSE)
