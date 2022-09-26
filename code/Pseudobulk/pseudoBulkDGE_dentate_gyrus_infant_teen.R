###########################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked dentate gyrus infant vs teen
# Anthony Ramnauth, Sept 26 2022
###########################################################

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
infant_teen_spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Adult")]
infant_teen_spe_pseudo <- infant_teen_spe_pseudo[, !infant_teen_spe_pseudo$age_bin %in% c("Elderly")]

# Format spe object for DE models
colData(infant_teen_spe_pseudo) <- colData(infant_teen_spe_pseudo)[, sort(c(
    "age",
    "age_bin",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells",
    "dentate_gyrus"
))]

colData(infant_teen_spe_pseudo)

colData(infant_teen_spe_pseudo)$ncells <- as.numeric(colData(infant_teen_spe_pseudo)$ncells)
colData(infant_teen_spe_pseudo)$race <- as.factor(colData(infant_teen_spe_pseudo)$race)
colData(infant_teen_spe_pseudo)$sample_id <- as.factor(colData(infant_teen_spe_pseudo)$sample_id)
colData(infant_teen_spe_pseudo)$sex <- as.factor(colData(infant_teen_spe_pseudo)$sex)
colData(infant_teen_spe_pseudo)$dentate_gyrus <- as.factor(colData(infant_teen_spe_pseudo)$dentate_gyrus)

colData(infant_teen_spe_pseudo)

# Drop things we don't need
spatialCoords(infant_teen_spe_pseudo) <- NULL
imgData(infant_teen_spe_pseudo) <- NULL

infant_teen_spe_pseudo$enrichment_infant <- 0
infant_teen_spe_pseudo$enrichment_infant[infant_teen_spe_pseudo$age_bin == "Infant"] <- 1

infant_teen_spe_pseudo$enrichment_teen <- 0
infant_teen_spe_pseudo$enrichment_teen[infant_teen_spe_pseudo$age_bin == "Teen"] <- 1


model_formula <- ~enrichment_infant
m <- model.matrix(model_formula, data = colData(infant_teen_spe_pseudo))

# Use pseudoBulkDGE to quickly perform age_bin DE analysis for BayesSpace labels
infant_de_results <- pseudoBulkDGE(
    infant_teen_spe_pseudo,
    col.data = colData(infant_teen_spe_pseudo),
    label = infant_teen_spe_pseudo$dentate_gyrus,
    design = ~enrichment_infant,
    coef = "enrichment_infant",
    row.data = rowData(infant_teen_spe_pseudo),
)

teen_de_results <- pseudoBulkDGE(
    infant_teen_spe_pseudo,
    col.data = colData(infant_teen_spe_pseudo),
    label = infant_teen_spe_pseudo$dentate_gyrus,
    design = ~enrichment_teen,
    coef = "enrichment_teen",
    row.data = rowData(infant_teen_spe_pseudo),
)

# Save modeling results
save(infant_de_results, teen_de_results,
    file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_dentategyrus_infant_teen_results.Rdata")
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
    "pseudoBulkDGE_DE_volcano_InfantvsTeen_dentategyrus.pdf"),
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
    subtitle = "Infant vs. Teen",
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
    subtitle = "Infant vs. Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

###########################################
# Volcano Plots of results for teen age_bin
###########################################

dentate0teen <- data.frame(
    gene_name = teen_de_results[[1]]$gene_name,
    logFC = teen_de_results[[1]]$logFC,
    FDR = teen_de_results[[1]]$FDR
)

dentate1teen <- data.frame(
    gene_name = teen_de_results[[2]]$gene_name,
    logFC = teen_de_results[[2]]$logFC,
    FDR = teen_de_results[[2]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE",
    "pseudoBulkDGE_DE_volcano_TeenvsInfant_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0teen,
    lab = dentate0teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 2,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Teen vs. Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1teen,
    lab = dentate1teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 2,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
        'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Teen vs. Infant",
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
