############################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked dentate gyrus infant vs adult
# Anthony Ramnauth, Sept 21 2022
############################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(sessioninfo)
    library(SingleCellExperiment)
    library(rafalib)
    library(limma)
    library(edgeR)
    library(scran)
    library(EnhancedVolcano)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Add colData() for Dentate Gyrus
spe_pseudo$dentate_gyrus <- 0
spe_pseudo$dentate_gyrus[spe_pseudo$BayesSpace == "1" |
        spe_pseudo$BayesSpace == "2"|
        spe_pseudo$BayesSpace == "4"|
        spe_pseudo$BayesSpace == "8"] <- 1

# Leave only infant vs. another age group
infant_adult_spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Teen")]
infant_adult_spe_pseudo <- infant_adult_spe_pseudo[, !infant_adult_spe_pseudo$age_bin %in% c("Elderly")]

# Format spe object for DE models
colData(infant_adult_spe_pseudo) <- colData(infant_adult_spe_pseudo)[, sort(c(
    "age",
    "age_bin",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells",
    "dentate_gyrus"
))]

colData(infant_adult_spe_pseudo)

colData(infant_adult_spe_pseudo)$ncells <- as.numeric(colData(infant_adult_spe_pseudo)$ncells)
colData(infant_adult_spe_pseudo)$race <- as.factor(colData(infant_adult_spe_pseudo)$race)
colData(infant_adult_spe_pseudo)$sample_id <- as.factor(colData(infant_adult_spe_pseudo)$sample_id)
colData(infant_adult_spe_pseudo)$sex <- as.factor(colData(infant_adult_spe_pseudo)$sex)
colData(infant_adult_spe_pseudo)$dentate_gyrus <- as.factor(colData(infant_adult_spe_pseudo)$dentate_gyrus)

colData(infant_adult_spe_pseudo)

# Drop things we don't need
spatialCoords(infant_adult_spe_pseudo) <- NULL
imgData(infant_adult_spe_pseudo) <- NULL

infant_adult_spe_pseudo$enrichment_infant <- 0
infant_adult_spe_pseudo$enrichment_infant[infant_adult_spe_pseudo$age_bin == "Infant"] <- 1

infant_adult_spe_pseudo$enrichment_adult <- 0
infant_adult_spe_pseudo$enrichment_adult[infant_adult_spe_pseudo$age_bin == "Adult"] <- 1


model_formula <- ~enrichment_infant
m <- model.matrix(model_formula, data = colData(infant_adult_spe_pseudo))

# Use pseudoBulkDGE to quickly perform age_bin DE analysis for BayesSpace labels
infant_de_results <- pseudoBulkDGE(
    infant_adult_spe_pseudo,
    col.data = colData(infant_adult_spe_pseudo),
    label = infant_adult_spe_pseudo$dentate_gyrus,
    design = ~enrichment_infant,
    coef = "enrichment_infant",
    row.data = rowData(infant_adult_spe_pseudo),
)

adult_de_results <- pseudoBulkDGE(
    infant_adult_spe_pseudo,
    col.data = colData(infant_adult_spe_pseudo),
    label = infant_adult_spe_pseudo$dentate_gyrus,
    design = ~enrichment_adult,
    coef = "enrichment_adult",
    row.data = rowData(infant_adult_spe_pseudo),
)

# Save modeling results
save(infant_de_results, adult_de_results,
    file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_dentategyrus_infant_adult_results.Rdata")
)

#############################################
# Volcano Plots of results for infant age_bin
#############################################

dentate0infant <- data.frame(
    gene_name = infant_de_results[[1]]$gene_name,
    logFC = infant_de_results[[1]]$logFC,
    FDR = infant_de_results[[1]]$FDR
)

dentate1infant <- data.frame(
    gene_name = infant_de_results[[2]]$gene_name,
    logFC = infant_de_results[[2]]$logFC,
    FDR = infant_de_results[[2]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE",
    "pseudoBulkDGE_DE_volcano_InfantvsAdult_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0infant,
    lab = dentate0infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 2,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1infant,
    lab = dentate1infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 2,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

############################################
# Volcano Plots of results for adult age_bin
############################################

dentate0adult <- data.frame(
    gene_name = adult_de_results[[1]]$gene_name,
    logFC = adult_de_results[[1]]$logFC,
    FDR = adult_de_results[[1]]$FDR
)

dentate1adult <- data.frame(
    gene_name = adult_de_results[[2]]$gene_name,
    logFC = adult_de_results[[2]]$logFC,
    FDR = adult_de_results[[2]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE",
    "pseudoBulkDGE_DE_volcano_AdultvsInfant_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0adult,
    lab = dentate0adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 2,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Adult vs. Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1adult,
    lab = dentate1adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 2,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
        'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Adult vs. Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
