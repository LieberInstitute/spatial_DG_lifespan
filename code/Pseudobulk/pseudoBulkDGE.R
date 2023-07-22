##############################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked across age_bin
# Anthony Ramnauth, June 13 2023
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
    library(dplyr)
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
colData(spe_pseudo)$BayesSpace <- as.factor(colData(spe_pseudo)$BayesSpace)


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

# Use pseudoBulkDGE to quickly perform age_bin DE analysis for BayesSpace labels
infant_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_infant,
    coef = "enrichment_infant",
    row.data = rowData(spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
)

teen_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_teen,
    coef = "enrichment_teen",
    row.data = rowData(spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
)

adult_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_adult,
    coef = "enrichment_adult",
    row.data = rowData(spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
)

elderly_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$BayesSpace,
    design = ~enrichment_elderly,
    coef = "enrichment_elderly",
    row.data = rowData(spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE)

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
    adj.P.Val = infant_de_results[[1]]$adj.P.Val
)

bayes10infant <- data.frame(
    gene_name = infant_de_results[[2]]$gene_name,
    logFC = infant_de_results[[2]]$logFC,
    adj.P.Val = infant_de_results[[2]]$adj.P.Val
)

bayes2infant <- data.frame(
    gene_name = infant_de_results[[3]]$gene_name,
    logFC = infant_de_results[[3]]$logFC,
    adj.P.Val = infant_de_results[[3]]$adj.P.Val
)

bayes4infant <- data.frame(
    gene_name = infant_de_results[[4]]$gene_name,
    logFC = infant_de_results[[4]]$logFC,
    adj.P.Val = infant_de_results[[4]]$adj.P.Val
)

bayes5infant <- data.frame(
    gene_name = infant_de_results[[5]]$gene_name,
    logFC = infant_de_results[[5]]$logFC,
    adj.P.Val = infant_de_results[[5]]$adj.P.Val
)

bayes6infant <- data.frame(
    gene_name = infant_de_results[[6]]$gene_name,
    logFC = infant_de_results[[6]]$logFC,
    adj.P.Val = infant_de_results[[6]]$adj.P.Val
)

bayes7infant <- data.frame(
    gene_name = infant_de_results[[7]]$gene_name,
    logFC = infant_de_results[[7]]$logFC,
    adj.P.Val = infant_de_results[[7]]$adj.P.Val
)

bayes8infant <- data.frame(
    gene_name = infant_de_results[[8]]$gene_name,
    logFC = infant_de_results[[8]]$logFC,
    adj.P.Val = infant_de_results[[8]]$adj.P.Val
)

bayes9infant <- data.frame(
    gene_name = infant_de_results[[9]]$gene_name,
    logFC = infant_de_results[[9]]$logFC,
    adj.P.Val = infant_de_results[[9]]$adj.P.Val
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Infant.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1infant,
    lab = bayes1infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes2infant,
    lab = bayes2infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes4infant,
    lab = bayes4infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes5infant,
    lab = bayes5infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes6infant,
    lab = bayes6infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes7infant,
    lab = bayes7infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes8infant,
    lab = bayes8infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes9infant,
    lab = bayes9infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 9",
    subtitle = "Infant vs. non-Infant"
    )

EnhancedVolcano(bayes10infant,
    lab = bayes10infant$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 10",
    subtitle = "Infant vs. non-Infant"
    )


dev.off()

###########################################
# Volcano Plots of results for teen age_bin
###########################################

bayes1teen <- data.frame(
    gene_name = teen_de_results[[1]]$gene_name,
    logFC = teen_de_results[[1]]$logFC,
    adj.P.Val = teen_de_results[[1]]$adj.P.Val
)

bayes10teen <- data.frame(
    gene_name = teen_de_results[[2]]$gene_name,
    logFC = teen_de_results[[2]]$logFC,
    adj.P.Val = teen_de_results[[2]]$adj.P.Val
)

bayes2teen <- data.frame(
    gene_name = teen_de_results[[3]]$gene_name,
    logFC = teen_de_results[[3]]$logFC,
    adj.P.Val = teen_de_results[[3]]$adj.P.Val
)

bayes4teen <- data.frame(
    gene_name = teen_de_results[[4]]$gene_name,
    logFC = teen_de_results[[4]]$logFC,
    adj.P.Val = teen_de_results[[4]]$adj.P.Val
)

bayes5teen <- data.frame(
    gene_name = teen_de_results[[5]]$gene_name,
    logFC = teen_de_results[[5]]$logFC,
    adj.P.Val = teen_de_results[[5]]$adj.P.Val
)

bayes6teen <- data.frame(
    gene_name = teen_de_results[[6]]$gene_name,
    logFC = teen_de_results[[6]]$logFC,
    adj.P.Val = teen_de_results[[6]]$adj.P.Val
)

bayes7teen <- data.frame(
    gene_name = teen_de_results[[7]]$gene_name,
    logFC = teen_de_results[[7]]$logFC,
    adj.P.Val = teen_de_results[[7]]$adj.P.Val
)

bayes8teen <- data.frame(
    gene_name = teen_de_results[[8]]$gene_name,
    logFC = teen_de_results[[8]]$logFC,
    adj.P.Val = teen_de_results[[8]]$adj.P.Val
)

bayes9teen <- data.frame(
    gene_name = teen_de_results[[9]]$gene_name,
    logFC = teen_de_results[[9]]$logFC,
    adj.P.Val = teen_de_results[[9]]$adj.P.Val
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Teen.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1teen,
    lab = bayes1teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes2teen,
    lab = bayes2teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
        'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes4teen,
    lab = bayes4teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes5teen,
    lab = bayes5teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes6teen,
    lab = bayes6teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes7teen,
    lab = bayes7teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes8teen,
    lab = bayes8teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes9teen,
    lab = bayes9teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 9",
    subtitle = "Teen vs. non-Teen"
    )

EnhancedVolcano(bayes10teen,
    lab = bayes10teen$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 10",
    subtitle = "Teen vs. non-Teen"
    )

dev.off()

############################################
# Volcano Plots of results for adult age_bin
############################################

bayes1adult <- data.frame(
    gene_name = adult_de_results[[1]]$gene_name,
    logFC = adult_de_results[[1]]$logFC,
    adj.P.Val = adult_de_results[[1]]$adj.P.Val
)

bayes10adult <- data.frame(
    gene_name = adult_de_results[[2]]$gene_name,
    logFC = adult_de_results[[2]]$logFC,
    adj.P.Val = adult_de_results[[2]]$adj.P.Val
)

bayes2adult <- data.frame(
    gene_name = adult_de_results[[3]]$gene_name,
    logFC = adult_de_results[[3]]$logFC,
    adj.P.Val = adult_de_results[[3]]$adj.P.Val
)

bayes4adult <- data.frame(
    gene_name = adult_de_results[[4]]$gene_name,
    logFC = adult_de_results[[4]]$logFC,
    adj.P.Val = adult_de_results[[4]]$adj.P.Val
)

bayes5adult <- data.frame(
    gene_name = adult_de_results[[5]]$gene_name,
    logFC = adult_de_results[[5]]$logFC,
    adj.P.Val = adult_de_results[[5]]$adj.P.Val
)

bayes6adult <- data.frame(
    gene_name = adult_de_results[[6]]$gene_name,
    logFC = adult_de_results[[6]]$logFC,
    adj.P.Val = adult_de_results[[6]]$adj.P.Val
)

bayes7adult <- data.frame(
    gene_name = adult_de_results[[7]]$gene_name,
    logFC = adult_de_results[[7]]$logFC,
    adj.P.Val = adult_de_results[[7]]$adj.P.Val
)

bayes8adult <- data.frame(
    gene_name = adult_de_results[[8]]$gene_name,
    logFC = adult_de_results[[8]]$logFC,
    adj.P.Val = adult_de_results[[8]]$adj.P.Val
)

bayes9adult <- data.frame(
    gene_name = adult_de_results[[9]]$gene_name,
    logFC = adult_de_results[[9]]$logFC,
    adj.P.Val = adult_de_results[[9]]$adj.P.Val
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Adult.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1adult,
    lab = bayes1adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 1",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes2adult,
    lab = bayes2adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes4adult,
    lab = bayes4adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes5adult,
    lab = bayes5adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes6adult,
    lab = bayes6adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes7adult,
    lab = bayes7adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes8adult,
    lab = bayes8adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes9adult,
    lab = bayes9adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 9",
    subtitle = "Adult vs. non-Adult"
    )

EnhancedVolcano(bayes10adult,
    lab = bayes10adult$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 10",
    subtitle = "Adult vs. non-Adult"
    )

dev.off()

##############################################
# Volcano Plots of results for elderly age_bin
##############################################

bayes1elderly <- data.frame(
    gene_name = elderly_de_results[[1]]$gene_name,
    logFC = elderly_de_results[[1]]$logFC,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val
)

bayes10elderly <- data.frame(
    gene_name = elderly_de_results[[2]]$gene_name,
    logFC = elderly_de_results[[2]]$logFC,
    adj.P.Val = elderly_de_results[[2]]$adj.P.Val
)

bayes2elderly <- data.frame(
    gene_name = elderly_de_results[[3]]$gene_name,
    logFC = elderly_de_results[[3]]$logFC,
    adj.P.Val = elderly_de_results[[3]]$adj.P.Val
)

bayes4elderly <- data.frame(
    gene_name = elderly_de_results[[4]]$gene_name,
    logFC = elderly_de_results[[4]]$logFC,
    adj.P.Val = elderly_de_results[[4]]$adj.P.Val
)

bayes5elderly <- data.frame(
    gene_name = elderly_de_results[[5]]$gene_name,
    logFC = elderly_de_results[[5]]$logFC,
    adj.P.Val = elderly_de_results[[5]]$adj.P.Val
)

bayes6elderly <- data.frame(
    gene_name = elderly_de_results[[6]]$gene_name,
    logFC = elderly_de_results[[6]]$logFC,
    adj.P.Val = elderly_de_results[[6]]$adj.P.Val
)

bayes7elderly <- data.frame(
    gene_name = elderly_de_results[[7]]$gene_name,
    logFC = elderly_de_results[[7]]$logFC,
    adj.P.Val = elderly_de_results[[7]]$adj.P.Val
)

bayes8elderly <- data.frame(
    gene_name = elderly_de_results[[8]]$gene_name,
    logFC = elderly_de_results[[8]]$logFC,
    adj.P.Val = elderly_de_results[[8]]$adj.P.Val
)

bayes9elderly <- data.frame(
    gene_name = elderly_de_results[[9]]$gene_name,
    logFC = elderly_de_results[[9]]$logFC,
    adj.P.Val = elderly_de_results[[9]]$adj.P.Val
)

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Elderly.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(bayes1elderly,
    lab = bayes1elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = c("HLA-DQA1", "HLA-DMA", "HLA-DMB", "HLA-DPB1", "CD14", "C1QC",
        "C1QB", "CHI3L1", "CCL2", "C1QA", "CHI3L2", "HCK", "CD68", "HAMP"),
    drawConnectors = TRUE,
    arrowheads = FALSE,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 1 ~ SLM",
    subtitle = "Elderly vs. non-Elderly"
    ) +
    xlim(c(-5, 5)) +
    ylim(c(0, 5))

EnhancedVolcano(bayes2elderly,
    lab = bayes2elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 2",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes4elderly,
    lab = bayes4elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 4",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes5elderly,
    lab = bayes5elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 5",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes6elderly,
    lab = bayes6elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 6",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes7elderly,
    lab = bayes7elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 7",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes8elderly,
    lab = bayes8elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 8",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes9elderly,
    lab = bayes9elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 9",
    subtitle = "Elderly vs. non-Elderly"
    )

EnhancedVolcano(bayes10elderly,
    lab = bayes10elderly$gene_name,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "BayesSpace cluster 10",
    subtitle = "Elderly vs. non-Elderly"
    )
dev.off()

######################################
# Write csv files for each DE analysis
######################################

# directory to save whole tissue results
dir_outputs <- here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results")

infant_bayes1 <- data.frame(
    gene_id = infant_de_results[[1]]$gene_id,
    gene_name = infant_de_results[[1]]$gene_name,
    gene_type = infant_de_results[[1]]$gene_type,
    pvalue = infant_de_results[[1]]$P.Value,
    adj.P.Val = infant_de_results[[1]]$adj.P.Val,
    logFC = infant_de_results[[1]]$logFC
)

infant_bayes1 <- infant_bayes1 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out1 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace1_DE")

# Export summary as .csv file
write.csv(infant_bayes1, fn_out1, row.names = FALSE)

infant_bayes10 <- data.frame(
    gene_id = infant_de_results[[2]]$gene_id,
    gene_name = infant_de_results[[2]]$gene_name,
    gene_type = infant_de_results[[2]]$gene_type,
    pvalue = infant_de_results[[2]]$P.Value,
    adj.P.Val = infant_de_results[[2]]$adj.P.Val,
    logFC = infant_de_results[[2]]$logFC
)

infant_bayes10 <- infant_bayes10 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out10 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace10_DE")

# Export summary as .csv file
write.csv(infant_bayes10, fn_out10, row.names = FALSE)

infant_bayes2 <- data.frame(
    gene_id = infant_de_results[[3]]$gene_id,
    gene_name = infant_de_results[[3]]$gene_name,
    gene_type = infant_de_results[[3]]$gene_type,
    pvalue = infant_de_results[[3]]$P.Value,
    adj.P.Val = infant_de_results[[3]]$adj.P.Val,
    logFC = infant_de_results[[3]]$logFC
)

infant_bayes2 <- infant_bayes2 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out2 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace2_DE")

# Export summary as .csv file
write.csv(infant_bayes2, fn_out2, row.names = FALSE)

infant_bayes4 <- data.frame(
    gene_id = infant_de_results[[4]]$gene_id,
    gene_name = infant_de_results[[4]]$gene_name,
    gene_type = infant_de_results[[4]]$gene_type,
    pvalue = infant_de_results[[4]]$P.Value,
    adj.P.Val = infant_de_results[[4]]$adj.P.Val,
    logFC = infant_de_results[[4]]$logFC
)

infant_bayes4 <- infant_bayes4 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out4 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace4_DE")

# Export summary as .csv file
write.csv(infant_bayes4, fn_out4, row.names = FALSE)

infant_bayes5 <- data.frame(
    gene_id = infant_de_results[[5]]$gene_id,
    gene_name = infant_de_results[[5]]$gene_name,
    gene_type = infant_de_results[[5]]$gene_type,
    pvalue = infant_de_results[[5]]$P.Value,
    adj.P.Val = infant_de_results[[5]]$adj.P.Val,
    logFC = infant_de_results[[5]]$logFC
)

infant_bayes5 <- infant_bayes5 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out5 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace5_DE")

# Export summary as .csv file
write.csv(infant_bayes5, fn_out5, row.names = FALSE)

infant_bayes6 <- data.frame(
    gene_id = infant_de_results[[6]]$gene_id,
    gene_name = infant_de_results[[6]]$gene_name,
    gene_type = infant_de_results[[6]]$gene_type,
    pvalue = infant_de_results[[6]]$P.Value,
    adj.P.Val = infant_de_results[[6]]$adj.P.Val,
    logFC = infant_de_results[[6]]$logFC
)

infant_bayes6 <- infant_bayes6 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out6 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace6_DE")

# Export summary as .csv file
write.csv(infant_bayes6, fn_out6, row.names = FALSE)

infant_bayes7 <- data.frame(
    gene_id = infant_de_results[[7]]$gene_id,
    gene_name = infant_de_results[[7]]$gene_name,
    gene_type = infant_de_results[[7]]$gene_type,
    pvalue = infant_de_results[[7]]$P.Value,
    adj.P.Val = infant_de_results[[7]]$adj.P.Val,
    logFC = infant_de_results[[7]]$logFC
)

infant_bayes7 <- infant_bayes7 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out7 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace7_DE")

# Export summary as .csv file
write.csv(infant_bayes7, fn_out7, row.names = FALSE)

infant_bayes8 <- data.frame(
    gene_id = infant_de_results[[8]]$gene_id,
    gene_name = infant_de_results[[8]]$gene_name,
    gene_type = infant_de_results[[8]]$gene_type,
    pvalue = infant_de_results[[8]]$P.Value,
    adj.P.Val = infant_de_results[[8]]$adj.P.Val,
    logFC = infant_de_results[[8]]$logFC
)

infant_bayes8 <- infant_bayes8 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out8 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace8_DE")

# Export summary as .csv file
write.csv(infant_bayes8, fn_out8, row.names = FALSE)

infant_bayes9 <- data.frame(
    gene_id = infant_de_results[[9]]$gene_id,
    gene_name = infant_de_results[[9]]$gene_name,
    gene_type = infant_de_results[[9]]$gene_type,
    pvalue = infant_de_results[[9]]$P.Value,
    adj.P.Val = infant_de_results[[9]]$adj.P.Val,
    logFC = infant_de_results[[9]]$logFC
)

infant_bayes9 <- infant_bayes9 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out9 <- file.path(dir_outputs, "InfantvsNonInfant_BayesSpace9_DE")

# Export summary as .csv file
write.csv(infant_bayes9, fn_out9, row.names = FALSE)

teen_bayes1 <- data.frame(
    gene_id = teen_de_results[[1]]$gene_id,
    gene_name = teen_de_results[[1]]$gene_name,
    gene_type = teen_de_results[[1]]$gene_type,
    pvalue = teen_de_results[[1]]$P.Value,
    adj.P.Val = teen_de_results[[1]]$adj.P.Val,
    logFC = teen_de_results[[1]]$logFC
)

teen_bayes1 <- teen_bayes1 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out11 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace1_DE")

# Export summary as .csv file
write.csv(teen_bayes1, fn_out11, row.names = FALSE)

teen_bayes10 <- data.frame(
    gene_id = teen_de_results[[2]]$gene_id,
    gene_name = teen_de_results[[2]]$gene_name,
    gene_type = teen_de_results[[2]]$gene_type,
    pvalue = teen_de_results[[2]]$P.Value,
    adj.P.Val = teen_de_results[[2]]$adj.P.Val,
    logFC = teen_de_results[[2]]$logFC
)

teen_bayes10 <- teen_bayes10 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out210 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace10_DE")

# Export summary as .csv file
write.csv(teen_bayes10, fn_out210, row.names = FALSE)

teen_bayes2 <- data.frame(
    gene_id = teen_de_results[[3]]$gene_id,
    gene_name = teen_de_results[[3]]$gene_name,
    gene_type = teen_de_results[[3]]$gene_type,
    pvalue = teen_de_results[[3]]$P.Value,
    adj.P.Val = teen_de_results[[3]]$adj.P.Val,
    logFC = teen_de_results[[3]]$logFC
)

teen_bayes2 <- teen_bayes2 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out22 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace2_DE")

# Export summary as .csv file
write.csv(teen_bayes2, fn_out22, row.names = FALSE)

teen_bayes4 <- data.frame(
    gene_id = teen_de_results[[4]]$gene_id,
    gene_name = teen_de_results[[4]]$gene_name,
    gene_type = teen_de_results[[4]]$gene_type,
    pvalue = teen_de_results[[4]]$P.Value,
    adj.P.Val = teen_de_results[[4]]$adj.P.Val,
    logFC = teen_de_results[[4]]$logFC
)

teen_bayes4 <- teen_bayes4 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out44 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace4_DE")

# Export summary as .csv file
write.csv(teen_bayes4, fn_out44, row.names = FALSE)

teen_bayes5 <- data.frame(
    gene_id = teen_de_results[[5]]$gene_id,
    gene_name = teen_de_results[[5]]$gene_name,
    gene_type = teen_de_results[[5]]$gene_type,
    pvalue = teen_de_results[[5]]$P.Value,
    adj.P.Val = teen_de_results[[5]]$adj.P.Val,
    logFC = teen_de_results[[5]]$logFC
)

teen_bayes5 <- teen_bayes5 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out55 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace5_DE")

# Export summary as .csv file
write.csv(teen_bayes5, fn_out55, row.names = FALSE)

teen_bayes6 <- data.frame(
    gene_id = teen_de_results[[6]]$gene_id,
    gene_name = teen_de_results[[6]]$gene_name,
    gene_type = teen_de_results[[6]]$gene_type,
    pvalue = teen_de_results[[6]]$P.Value,
    adj.P.Val = teen_de_results[[6]]$adj.P.Val,
    logFC = teen_de_results[[6]]$logFC
)

teen_bayes6 <- teen_bayes6 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out66 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace6_DE")

# Export summary as .csv file
write.csv(teen_bayes6, fn_out66, row.names = FALSE)

teen_bayes7 <- data.frame(
    gene_id = teen_de_results[[7]]$gene_id,
    gene_name = teen_de_results[[7]]$gene_name,
    gene_type = teen_de_results[[7]]$gene_type,
    pvalue = teen_de_results[[7]]$P.Value,
    adj.P.Val = teen_de_results[[7]]$adj.P.Val,
    logFC = teen_de_results[[7]]$logFC
)

teen_bayes7 <- teen_bayes7 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out77 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace7_DE")

# Export summary as .csv file
write.csv(teen_bayes7, fn_out77, row.names = FALSE)

teen_bayes8 <- data.frame(
    gene_id = teen_de_results[[8]]$gene_id,
    gene_name = teen_de_results[[8]]$gene_name,
    gene_type = teen_de_results[[8]]$gene_type,
    pvalue = teen_de_results[[8]]$P.Value,
    adj.P.Val = teen_de_results[[8]]$adj.P.Val,
    logFC = teen_de_results[[8]]$logFC
)

teen_bayes8 <- teen_bayes8 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out88 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace8_DE")

# Export summary as .csv file
write.csv(teen_bayes8, fn_out88, row.names = FALSE)

teen_bayes9 <- data.frame(
    gene_id = teen_de_results[[9]]$gene_id,
    gene_name = teen_de_results[[9]]$gene_name,
    gene_type = teen_de_results[[9]]$gene_type,
    pvalue = teen_de_results[[9]]$P.Value,
    adj.P.Val = teen_de_results[[9]]$adj.P.Val,
    logFC = teen_de_results[[9]]$logFC
)

teen_bayes9 <- teen_bayes9 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out99 <- file.path(dir_outputs, "TeenvsNonTeen_BayesSpace9_DE")

# Export summary as .csv file
write.csv(teen_bayes9, fn_out99, row.names = FALSE)

adult_bayes1 <- data.frame(
    gene_id = adult_de_results[[1]]$gene_id,
    gene_name = adult_de_results[[1]]$gene_name,
    gene_type = adult_de_results[[1]]$gene_type,
    pvalue = adult_de_results[[1]]$P.Value,
    adj.P.Val = adult_de_results[[1]]$adj.P.Val,
    logFC = adult_de_results[[1]]$logFC
)

adult_bayes1 <- adult_bayes1 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out111 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace1_DE")

# Export summary as .csv file
write.csv(adult_bayes1, fn_out111, row.names = FALSE)

adult_bayes10 <- data.frame(
    gene_id = adult_de_results[[2]]$gene_id,
    gene_name = adult_de_results[[2]]$gene_name,
    gene_type = adult_de_results[[2]]$gene_type,
    pvalue = adult_de_results[[2]]$P.Value,
    adj.P.Val = adult_de_results[[2]]$adj.P.Val,
    logFC = adult_de_results[[2]]$logFC
)

adult_bayes10 <- adult_bayes10 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out310 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace10_DE")

# Export summary as .csv file
write.csv(adult_bayes10, fn_out310, row.names = FALSE)

adult_bayes2 <- data.frame(
    gene_id = adult_de_results[[3]]$gene_id,
    gene_name = adult_de_results[[3]]$gene_name,
    gene_type = adult_de_results[[3]]$gene_type,
    pvalue = adult_de_results[[3]]$P.Value,
    adj.P.Val = adult_de_results[[3]]$adj.P.Val,
    logFC = adult_de_results[[3]]$logFC
)

adult_bayes2 <- adult_bayes2 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out222 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace2_DE")

# Export summary as .csv file
write.csv(adult_bayes2, fn_out222, row.names = FALSE)

adult_bayes4 <- data.frame(
    gene_id = adult_de_results[[4]]$gene_id,
    gene_name = adult_de_results[[4]]$gene_name,
    gene_type = adult_de_results[[4]]$gene_type,
    pvalue = adult_de_results[[4]]$P.Value,
    adj.P.Val = adult_de_results[[4]]$adj.P.Val,
    logFC = adult_de_results[[4]]$logFC
)

adult_bayes4 <- adult_bayes4 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out444 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace4_DE")

# Export summary as .csv file
write.csv(adult_bayes4, fn_out444, row.names = FALSE)

adult_bayes5 <- data.frame(
    gene_id = adult_de_results[[5]]$gene_id,
    gene_name = adult_de_results[[5]]$gene_name,
    gene_type = adult_de_results[[5]]$gene_type,
    pvalue = adult_de_results[[5]]$P.Value,
    adj.P.Val = adult_de_results[[5]]$adj.P.Val,
    logFC = adult_de_results[[5]]$logFC
)

adult_bayes5 <- adult_bayes5 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out555 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace5_DE")

# Export summary as .csv file
write.csv(adult_bayes5, fn_out555, row.names = FALSE)

adult_bayes6 <- data.frame(
    gene_id = adult_de_results[[6]]$gene_id,
    gene_name = adult_de_results[[6]]$gene_name,
    gene_type = adult_de_results[[6]]$gene_type,
    pvalue = adult_de_results[[6]]$P.Value,
    adj.P.Val = adult_de_results[[6]]$adj.P.Val,
    logFC = adult_de_results[[6]]$logFC
)

adult_bayes6 <- adult_bayes6 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out666 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace6_DE")

# Export summary as .csv file
write.csv(adult_bayes6, fn_out666, row.names = FALSE)

adult_bayes7 <- data.frame(
    gene_id = adult_de_results[[7]]$gene_id,
    gene_name = adult_de_results[[7]]$gene_name,
    gene_type = adult_de_results[[7]]$gene_type,
    pvalue = adult_de_results[[7]]$P.Value,
    adj.P.Val = adult_de_results[[7]]$adj.P.Val,
    logFC = adult_de_results[[7]]$logFC
)

adult_bayes7 <- adult_bayes7 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out777 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace7_DE")

# Export summary as .csv file
write.csv(adult_bayes7, fn_out777, row.names = FALSE)

adult_bayes8 <- data.frame(
    gene_id = adult_de_results[[8]]$gene_id,
    gene_name = adult_de_results[[8]]$gene_name,
    gene_type = adult_de_results[[8]]$gene_type,
    pvalue = adult_de_results[[8]]$P.Value,
    adj.P.Val = adult_de_results[[8]]$adj.P.Val,
    logFC = adult_de_results[[8]]$logFC
)

adult_bayes8 <- adult_bayes8 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out888 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace8_DE")

# Export summary as .csv file
write.csv(adult_bayes8, fn_out888, row.names = FALSE)

adult_bayes9 <- data.frame(
    gene_id = adult_de_results[[9]]$gene_id,
    gene_name = adult_de_results[[9]]$gene_name,
    gene_type = adult_de_results[[9]]$gene_type,
    pvalue = adult_de_results[[9]]$P.Value,
    adj.P.Val = adult_de_results[[9]]$adj.P.Val,
    logFC = adult_de_results[[9]]$logFC
)

adult_bayes9 <- adult_bayes9 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out999 <- file.path(dir_outputs, "AdultvsNonAdult_BayesSpace9_DE")

# Export summary as .csv file
write.csv(adult_bayes9, fn_out999, row.names = FALSE)

elderly_bayes1 <- data.frame(
    gene_id = elderly_de_results[[1]]$gene_id,
    gene_name = elderly_de_results[[1]]$gene_name,
    gene_type = elderly_de_results[[1]]$gene_type,
    pvalue = elderly_de_results[[1]]$P.Value,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val,
    logFC = elderly_de_results[[1]]$logFC
)

elderly_bayes1 <- elderly_bayes1 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out1111 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace1_DE")

# Export summary as .csv file
write.csv(elderly_bayes1, fn_out1111, row.names = FALSE)

elderly_bayes10 <- data.frame(
    gene_id = elderly_de_results[[2]]$gene_id,
    gene_name = elderly_de_results[[2]]$gene_name,
    gene_type = elderly_de_results[[2]]$gene_type,
    pvalue = elderly_de_results[[2]]$P.Value,
    adj.P.Val = elderly_de_results[[2]]$adj.P.Val,
    logFC = elderly_de_results[[2]]$logFC
)

elderly_bayes10 <- elderly_bayes10 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out410 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace10_DE")

# Export summary as .csv file
write.csv(elderly_bayes10, fn_out410, row.names = FALSE)

elderly_bayes2 <- data.frame(
    gene_id = elderly_de_results[[3]]$gene_id,
    gene_name = elderly_de_results[[3]]$gene_name,
    gene_type = elderly_de_results[[3]]$gene_type,
    pvalue = elderly_de_results[[3]]$P.Value,
    adj.P.Val = elderly_de_results[[3]]$adj.P.Val,
    logFC = elderly_de_results[[3]]$logFC
)

elderly_bayes2 <- elderly_bayes2 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out2222 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace2_DE")

# Export summary as .csv file
write.csv(elderly_bayes2, fn_out2222, row.names = FALSE)

elderly_bayes4 <- data.frame(
    gene_id = elderly_de_results[[4]]$gene_id,
    gene_name = elderly_de_results[[4]]$gene_name,
    gene_type = elderly_de_results[[4]]$gene_type,
    pvalue = elderly_de_results[[4]]$P.Value,
    adj.P.Val = elderly_de_results[[4]]$adj.P.Val,
    logFC = elderly_de_results[[4]]$logFC
)

elderly_bayes4 <- elderly_bayes4 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out4444 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace4_DE")

# Export summary as .csv file
write.csv(elderly_bayes4, fn_out4444, row.names = FALSE)

elderly_bayes5 <- data.frame(
    gene_id = elderly_de_results[[5]]$gene_id,
    gene_name = elderly_de_results[[5]]$gene_name,
    gene_type = elderly_de_results[[5]]$gene_type,
    pvalue = elderly_de_results[[5]]$P.Value,
    adj.P.Val = elderly_de_results[[5]]$adj.P.Val,
    logFC = elderly_de_results[[5]]$logFC
)

elderly_bayes5 <- elderly_bayes5 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out5555 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace5_DE")

# Export summary as .csv file
write.csv(elderly_bayes5, fn_out5555, row.names = FALSE)

elderly_bayes6 <- data.frame(
    gene_id = elderly_de_results[[6]]$gene_id,
    gene_name = elderly_de_results[[6]]$gene_name,
    gene_type = elderly_de_results[[6]]$gene_type,
    pvalue = elderly_de_results[[6]]$P.Value,
    adj.P.Val = elderly_de_results[[6]]$adj.P.Val,
    logFC = elderly_de_results[[6]]$logFC
)

elderly_bayes6 <- elderly_bayes6 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out6666 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace6_DE")

# Export summary as .csv file
write.csv(elderly_bayes6, fn_out6666, row.names = FALSE)

elderly_bayes7 <- data.frame(
    gene_id = elderly_de_results[[7]]$gene_id,
    gene_name = elderly_de_results[[7]]$gene_name,
    gene_type = elderly_de_results[[7]]$gene_type,
    pvalue = elderly_de_results[[7]]$P.Value,
    adj.P.Val = elderly_de_results[[7]]$adj.P.Val,
    logFC = elderly_de_results[[7]]$logFC
)

elderly_bayes7 <- elderly_bayes7 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out7777 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace7_DE")

# Export summary as .csv file
write.csv(elderly_bayes7, fn_out7777, row.names = FALSE)

elderly_bayes8 <- data.frame(
    gene_id = elderly_de_results[[8]]$gene_id,
    gene_name = elderly_de_results[[8]]$gene_name,
    gene_type = elderly_de_results[[8]]$gene_type,
    pvalue = elderly_de_results[[8]]$P.Value,
    adj.P.Val = elderly_de_results[[8]]$adj.P.Val,
    logFC = elderly_de_results[[8]]$logFC
)

elderly_bayes8 <- elderly_bayes8 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out8888 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace8_DE")

# Export summary as .csv file
write.csv(elderly_bayes8, fn_out8888, row.names = FALSE)

elderly_bayes9 <- data.frame(
    gene_id = elderly_de_results[[9]]$gene_id,
    gene_name = elderly_de_results[[9]]$gene_name,
    gene_type = elderly_de_results[[9]]$gene_type,
    pvalue = elderly_de_results[[9]]$P.Value,
    adj.P.Val = elderly_de_results[[9]]$adj.P.Val,
    logFC = elderly_de_results[[9]]$logFC
)

elderly_bayes9 <- elderly_bayes9 %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out9999 <- file.path(dir_outputs, "ElderlyvsNonElderly_BayesSpace9_DE")

# Export summary as .csv file
write.csv(elderly_bayes9, fn_out9999, row.names = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
