###########################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked dentate gyrus across age_bin
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

# Format spe object for DE models
colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
    "age",
    "age_bin",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells",
    "dentate_gyrus"
))]

colData(spe_pseudo)

colData(spe_pseudo)$ncells <- as.numeric(colData(spe_pseudo)$ncells)
colData(spe_pseudo)$race <- as.factor(colData(spe_pseudo)$race)
colData(spe_pseudo)$sample_id <- as.factor(colData(spe_pseudo)$sample_id)
colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)
colData(spe_pseudo)$dentate_gyrus <- as.factor(colData(spe_pseudo)$dentate_gyrus)


colData(spe_pseudo)

# Drop things we don't need
spatialCoords(spe_pseudo) <- NULL
imgData(spe_pseudo) <- NULL

spe_pseudo$enrichment_infant <- 0
spe_pseudo$enrichment_infant[spe_pseudo$age_bin == "Infant"] <- 1

spe_pseudo$enrichment_teen <- 0
spe_pseudo$enrichment_teen[spe_pseudo$age_bin == "Teen"] <- 1

spe_pseudo$enrichment_adult <- 0
spe_pseudo$enrichment_adult[spe_pseudo$age_bin == "Adult"] <- 1

spe_pseudo$enrichment_elderly <- 0
spe_pseudo$enrichment_elderly[spe_pseudo$age_bin == "Elderly"] <- 1

model_formula <- ~enrichment_elderly
m <- model.matrix(model_formula, data = colData(spe_pseudo))

# Use pseudoBulkDGE to quickly perform age_bin DE analysis for BayesSpace labels
infant_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$dentate_gyrus,
    design = ~enrichment_infant,
    coef = "enrichment_infant",
    row.data = rowData(spe_pseudo),
)

teen_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$dentate_gyrus,
    design = ~enrichment_teen,
    coef = "enrichment_teen",
    row.data = rowData(spe_pseudo),
)

adult_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$dentate_gyrus,
    design = ~enrichment_adult,
    coef = "enrichment_adult",
    row.data = rowData(spe_pseudo),
)

elderly_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$dentate_gyrus,
    design = ~enrichment_elderly,
    coef = "enrichment_elderly",
    row.data = rowData(spe_pseudo),
)

# Save modeling results
save(infant_de_results, teen_de_results, adult_de_results, elderly_de_results,
    file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_dentategyrus_results.Rdata")
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

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Infant_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0infant,
    lab = dentate0infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1infant,
    lab = dentate1infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Infant vs. non-Infant",
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

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Teen_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0teen,
    lab = dentate0teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1teen,
    lab = dentate1teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
        'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Teen vs. non-Teen",
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

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Adult_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0adult,
    lab = dentate0adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1adult,
    lab = dentate1adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

##############################################
# Volcano Plots of results for elderly age_bin
##############################################

dentate0elderly <- data.frame(
    gene_name = elderly_de_results[[1]]$gene_name,
    logFC = elderly_de_results[[1]]$logFC,
    FDR = elderly_de_results[[1]]$FDR
)

dentate1elderly <- data.frame(
    gene_name = elderly_de_results[[2]]$gene_name,
    logFC = elderly_de_results[[2]]$logFC,
    FDR = elderly_de_results[[2]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Elderly_dentategyrus.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(dentate0elderly,
    lab = dentate0elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Non-Dentate Gyrus",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(dentate1elderly,
    lab = dentate1elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Elderly vs. non-Elderly",
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
