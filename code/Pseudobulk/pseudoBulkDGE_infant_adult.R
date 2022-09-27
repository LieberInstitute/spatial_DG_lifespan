##############################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked infant vs adult
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
    "ncells"
))]

colData(infant_adult_spe_pseudo)

colData(infant_adult_spe_pseudo)$ncells <- as.numeric(colData(infant_adult_spe_pseudo)$ncells)
colData(infant_adult_spe_pseudo)$race <- as.factor(colData(infant_adult_spe_pseudo)$race)
colData(infant_adult_spe_pseudo)$sample_id <- as.factor(colData(infant_adult_spe_pseudo)$sample_id)
colData(infant_adult_spe_pseudo)$sex <- as.factor(colData(infant_adult_spe_pseudo)$sex)

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
    label = infant_adult_spe_pseudo$BayesSpace,
    design = ~enrichment_infant,
    coef = "enrichment_infant",
    row.data = rowData(infant_adult_spe_pseudo),
)

adult_de_results <- pseudoBulkDGE(
    infant_adult_spe_pseudo,
    col.data = colData(infant_adult_spe_pseudo),
    label = infant_adult_spe_pseudo$BayesSpace,
    design = ~enrichment_adult,
    coef = "enrichment_adult",
    row.data = rowData(infant_adult_spe_pseudo),
)

# Save modeling results
save(infant_de_results, adult_de_results,
    file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_infant_adult_results.Rdata")
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

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_InfantvsAdult.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1infant,
    lab = bayes1infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes2infant,
    lab = bayes2infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes3infant,
    lab = bayes3infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 3",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes4infant,
    lab = bayes4infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes5infant,
    lab = bayes5infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes6infant,
    lab = bayes6infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes7infant,
    lab = bayes7infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

EnhancedVolcano(bayes8infant,
    lab = bayes8infant$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Infant vs. Adult",
    drawConnectors = TRUE,
    colConnectors = 'black'
    )

dev.off()

###########################################
# Volcano Plots of results for teen age_bin
###########################################

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

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_AdultvsInfant.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1adult,
    lab = bayes1adult$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    ylab = "-log10 FDR",
    pCutoff = 0.049,
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Adult vs. Infant",
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
    subtitle = "Adult vs. Infant",
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
    subtitle = "Adult vs. Infant",
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
    subtitle = "Adult vs. Infant",
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
    subtitle = "Adult vs. Infant",
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
    subtitle = "Adult vs. Infant",
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
    subtitle = "Adult vs. Infant",
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
