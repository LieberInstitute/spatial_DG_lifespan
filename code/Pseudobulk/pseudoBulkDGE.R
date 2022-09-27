##############################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked across age_bin
# Anthony Ramnauth, Sept 21 2022
##############################################

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

# Format spe object for DE models

colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
    "age",
    "age_bin",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells"
))]

colData(spe_pseudo)

colData(spe_pseudo)$ncells <- as.numeric(colData(spe_pseudo)$ncells)
colData(spe_pseudo)$race <- as.factor(colData(spe_pseudo)$race)
colData(spe_pseudo)$sample_id <- as.factor(colData(spe_pseudo)$sample_id)
colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)

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
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_infant,
    coef = "enrichment_infant",
    row.data = rowData(spe_pseudo),
)

teen_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_teen,
    coef = "enrichment_teen",
    row.data = rowData(spe_pseudo),
)

adult_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_adult,
    coef = "enrichment_adult",
    row.data = rowData(spe_pseudo),
)

elderly_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_elderly,
    coef = "enrichment_elderly",
    row.data = rowData(spe_pseudo),
)

# Save modeling results
save(infant_de_results, teen_de_results, adult_de_results, elderly_de_results,
    file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results.Rdata")
)

#############################################
# Volcano Plots of results for infant age_bin
#############################################

bayes1infant <- data.frame(
    gene_name = infant_de_results[[1]]$gene_name,
    logFC = infant_de_results[[1]]$logFC,
    FDR = infant_de_results[[1]]$FDR
)

bayes2infant <- data.frame(
    gene_name = infant_de_results[[2]]$gene_name,
    logFC = infant_de_results[[2]]$logFC,
    FDR = infant_de_results[[2]]$FDR
)

bayes3infant <- data.frame(
    gene_name = infant_de_results[[3]]$gene_name,
    logFC = infant_de_results[[3]]$logFC,
    FDR = infant_de_results[[3]]$FDR
)

bayes4infant <- data.frame(
    gene_name = infant_de_results[[4]]$gene_name,
    logFC = infant_de_results[[4]]$logFC,
    FDR = infant_de_results[[4]]$FDR
)

bayes5infant <- data.frame(
    gene_name = infant_de_results[[5]]$gene_name,
    logFC = infant_de_results[[5]]$logFC,
    FDR = infant_de_results[[5]]$FDR
)

bayes6infant <- data.frame(
    gene_name = infant_de_results[[6]]$gene_name,
    logFC = infant_de_results[[6]]$logFC,
    FDR = infant_de_results[[6]]$FDR
)

bayes7infant <- data.frame(
    gene_name = infant_de_results[[7]]$gene_name,
    logFC = infant_de_results[[7]]$logFC,
    FDR = infant_de_results[[7]]$FDR
)

bayes8infant <- data.frame(
    gene_name = infant_de_results[[8]]$gene_name,
    logFC = infant_de_results[[8]]$logFC,
    FDR = infant_de_results[[8]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Infant.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1infant,
    lab = bayes1infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes2infant,
    lab = bayes2infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes3infant,
    lab = bayes3infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 3",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes4infant,
    lab = bayes4infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes5infant,
    lab = bayes5infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes6infant,
    lab = bayes6infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes7infant,
    lab = bayes7infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes8infant,
    lab = bayes8infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Infant vs. non-Infant",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

###########################################
# Volcano Plots of results for teen age_bin
###########################################

bayes1teen <- data.frame(
    gene_name = teen_de_results[[1]]$gene_name,
    logFC = teen_de_results[[1]]$logFC,
    FDR = teen_de_results[[1]]$FDR
)

bayes2teen <- data.frame(
    gene_name = teen_de_results[[2]]$gene_name,
    logFC = teen_de_results[[2]]$logFC,
    FDR = teen_de_results[[2]]$FDR
)

bayes3teen <- data.frame(
    gene_name = teen_de_results[[3]]$gene_name,
    logFC = teen_de_results[[3]]$logFC,
    FDR = teen_de_results[[3]]$FDR
)

bayes4teen <- data.frame(
    gene_name = teen_de_results[[4]]$gene_name,
    logFC = teen_de_results[[4]]$logFC,
    FDR = teen_de_results[[4]]$FDR
)

bayes5teen <- data.frame(
    gene_name = teen_de_results[[5]]$gene_name,
    logFC = teen_de_results[[5]]$logFC,
    FDR = teen_de_results[[5]]$FDR
)

bayes6teen <- data.frame(
    gene_name = teen_de_results[[6]]$gene_name,
    logFC = teen_de_results[[6]]$logFC,
    FDR = teen_de_results[[6]]$FDR
)

bayes7teen <- data.frame(
    gene_name = teen_de_results[[7]]$gene_name,
    logFC = teen_de_results[[7]]$logFC,
    FDR = teen_de_results[[7]]$FDR
)

bayes8teen <- data.frame(
    gene_name = teen_de_results[[8]]$gene_name,
    logFC = teen_de_results[[8]]$logFC,
    FDR = teen_de_results[[8]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Teen.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1teen,
    lab = bayes1teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes2teen,
    lab = bayes2teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
        'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes3teen,
    lab = bayes3teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 3",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes4teen,
    lab = bayes4teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes5teen,
    lab = bayes5teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes6teen,
    lab = bayes6teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes7teen,
    lab = bayes7teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes8teen,
    lab = bayes8teen$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Teen vs. non-Teen",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

############################################
# Volcano Plots of results for adult age_bin
############################################

bayes1adult <- data.frame(
    gene_name = adult_de_results[[1]]$gene_name,
    logFC = adult_de_results[[1]]$logFC,
    FDR = adult_de_results[[1]]$FDR
)

bayes2adult <- data.frame(
    gene_name = adult_de_results[[2]]$gene_name,
    logFC = adult_de_results[[2]]$logFC,
    FDR = adult_de_results[[2]]$FDR
)

bayes3adult <- data.frame(
    gene_name = adult_de_results[[3]]$gene_name,
    logFC = adult_de_results[[3]]$logFC,
    FDR = adult_de_results[[3]]$FDR
)

bayes4adult <- data.frame(
    gene_name = adult_de_results[[4]]$gene_name,
    logFC = adult_de_results[[4]]$logFC,
    FDR = adult_de_results[[4]]$FDR
)

bayes5adult <- data.frame(
    gene_name = adult_de_results[[5]]$gene_name,
    logFC = adult_de_results[[5]]$logFC,
    FDR = adult_de_results[[5]]$FDR
)

bayes6adult <- data.frame(
    gene_name = adult_de_results[[6]]$gene_name,
    logFC = adult_de_results[[6]]$logFC,
    FDR = adult_de_results[[6]]$FDR
)

bayes7adult <- data.frame(
    gene_name = adult_de_results[[7]]$gene_name,
    logFC = adult_de_results[[7]]$logFC,
    FDR = adult_de_results[[7]]$FDR
)

bayes8adult <- data.frame(
    gene_name = adult_de_results[[8]]$gene_name,
    logFC = adult_de_results[[8]]$logFC,
    FDR = adult_de_results[[8]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Adult.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1adult,
    lab = bayes1adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes2adult,
    lab = bayes2adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes3adult,
    lab = bayes3adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 3",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes4adult,
    lab = bayes4adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes5adult,
    lab = bayes5adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes6adult,
    lab = bayes6adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes7adult,
    lab = bayes7adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes8adult,
    lab = bayes8adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Adult vs. non-Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

##############################################
# Volcano Plots of results for elderly age_bin
##############################################

bayes1elderly <- data.frame(
    gene_name = elderly_de_results[[1]]$gene_name,
    logFC = elderly_de_results[[1]]$logFC,
    FDR = elderly_de_results[[1]]$FDR
)

bayes2elderly <- data.frame(
    gene_name = elderly_de_results[[2]]$gene_name,
    logFC = elderly_de_results[[2]]$logFC,
    FDR = elderly_de_results[[2]]$FDR
)

bayes3elderly <- data.frame(
    gene_name = elderly_de_results[[3]]$gene_name,
    logFC = elderly_de_results[[3]]$logFC,
    FDR = elderly_de_results[[3]]$FDR
)

bayes4elderly <- data.frame(
    gene_name = elderly_de_results[[4]]$gene_name,
    logFC = elderly_de_results[[4]]$logFC,
    FDR = elderly_de_results[[4]]$FDR
)

bayes5elderly <- data.frame(
    gene_name = elderly_de_results[[5]]$gene_name,
    logFC = elderly_de_results[[5]]$logFC,
    FDR = elderly_de_results[[5]]$FDR
)

bayes6elderly <- data.frame(
    gene_name = elderly_de_results[[6]]$gene_name,
    logFC = elderly_de_results[[6]]$logFC,
    FDR = elderly_de_results[[6]]$FDR
)

bayes7elderly <- data.frame(
    gene_name = elderly_de_results[[7]]$gene_name,
    logFC = elderly_de_results[[7]]$logFC,
    FDR = elderly_de_results[[7]]$FDR
)

bayes8elderly <- data.frame(
    gene_name = elderly_de_results[[8]]$gene_name,
    logFC = elderly_de_results[[8]]$logFC,
    FDR = elderly_de_results[[8]]$FDR
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Elderly.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1elderly,
    lab = bayes1elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes2elderly,
    lab = bayes2elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes3elderly,
    lab = bayes3elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 3",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes4elderly,
    lab = bayes4elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes5elderly,
    lab = bayes5elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes6elderly,
    lab = bayes6elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes7elderly,
    lab = bayes7elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Elderly vs. non-Elderly",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes8elderly,
    lab = bayes8elderly$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
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
