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
    library(dplyr)
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

######################################
# Write csv files for each DE analysis
######################################

# directory to save whole tissue results
dir_outputs <- here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results")

infant_dg0 <- data.frame(
    gene_id = infant_de_results[[1]]$gene_id,
    gene_name = infant_de_results[[1]]$gene_name,
    gene_type = infant_de_results[[1]]$gene_type,
    pvalue = infant_de_results[[1]]$PValue,
    FDR = infant_de_results[[1]]$FDR,
    logFC = infant_de_results[[1]]$logFC
)

infant_dg0 <- infant_dg0 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out1 <- file.path(dir_outputs, "InfantvsNonInfant_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(infant_dg0, fn_out1, row.names = FALSE)

infant_dg1 <- data.frame(
    gene_id = infant_de_results[[2]]$gene_id,
    gene_name = infant_de_results[[2]]$gene_name,
    gene_type = infant_de_results[[2]]$gene_type,
    pvalue = infant_de_results[[2]]$PValue,
    FDR = infant_de_results[[2]]$FDR,
    logFC = infant_de_results[[2]]$logFC
)

infant_dg1 <- infant_dg1 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out2 <- file.path(dir_outputs, "InfantvsNonInfant_DentateGyrus_DE")

# Export summary as .csv file
write.csv(infant_dg1, fn_out2, row.names = FALSE)

teen_dg0 <- data.frame(
    gene_id = teen_de_results[[1]]$gene_id,
    gene_name = teen_de_results[[1]]$gene_name,
    gene_type = teen_de_results[[1]]$gene_type,
    pvalue = teen_de_results[[1]]$PValue,
    FDR = teen_de_results[[1]]$FDR,
    logFC = teen_de_results[[1]]$logFC
)

teen_dg0 <- teen_dg0 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out3 <- file.path(dir_outputs, "TeenvsNonTeen_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(teen_dg0, fn_out3, row.names = FALSE)

teen_dg1 <- data.frame(
    gene_id = teen_de_results[[2]]$gene_id,
    gene_name = teen_de_results[[2]]$gene_name,
    gene_type = teen_de_results[[2]]$gene_type,
    pvalue = teen_de_results[[2]]$PValue,
    FDR = teen_de_results[[2]]$FDR,
    logFC = teen_de_results[[2]]$logFC
)

teen_dg1 <- teen_dg1 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out4 <- file.path(dir_outputs, "TeenvsNonTeen_DentateGyrus_DE")

# Export summary as .csv file
write.csv(teen_dg1, fn_out4, row.names = FALSE)

adult_dg0 <- data.frame(
    gene_id = adult_de_results[[1]]$gene_id,
    gene_name = adult_de_results[[1]]$gene_name,
    gene_type = adult_de_results[[1]]$gene_type,
    pvalue = adult_de_results[[1]]$PValue,
    FDR = adult_de_results[[1]]$FDR,
    logFC = adult_de_results[[1]]$logFC
)

adult_dg0 <- adult_dg0 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out5 <- file.path(dir_outputs, "AdultvsNonAdult_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(adult_dg0, fn_out5, row.names = FALSE)

adult_dg1 <- data.frame(
    gene_id = adult_de_results[[2]]$gene_id,
    gene_name = adult_de_results[[2]]$gene_name,
    gene_type = adult_de_results[[2]]$gene_type,
    pvalue = adult_de_results[[2]]$PValue,
    FDR = adult_de_results[[2]]$FDR,
    logFC = adult_de_results[[2]]$logFC
)

adult_dg1 <- adult_dg1 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out6 <- file.path(dir_outputs, "AdultvsNonAdult_DentateGyrus_DE")

# Export summary as .csv file
write.csv(adult_dg1, fn_out6, row.names = FALSE)

elderly_dg0 <- data.frame(
    gene_id = elderly_de_results[[1]]$gene_id,
    gene_name = elderly_de_results[[1]]$gene_name,
    gene_type = elderly_de_results[[1]]$gene_type,
    pvalue = elderly_de_results[[1]]$PValue,
    FDR = elderly_de_results[[1]]$FDR,
    logFC = elderly_de_results[[1]]$logFC
)

elderly_dg0 <- elderly_dg0 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out7 <- file.path(dir_outputs, "ElderlyvsNonElderly_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(elderly_dg0, fn_out7, row.names = FALSE)

elderly_dg1 <- data.frame(
    gene_id = elderly_de_results[[2]]$gene_id,
    gene_name = elderly_de_results[[2]]$gene_name,
    gene_type = elderly_de_results[[2]]$gene_type,
    pvalue = elderly_de_results[[2]]$PValue,
    FDR = elderly_de_results[[2]]$FDR,
    logFC = elderly_de_results[[2]]$logFC
)

elderly_dg1 <- elderly_dg1 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out8 <- file.path(dir_outputs, "ElderlyvsNonElderly_DentateGyrus_DE")

# Export summary as .csv file
write.csv(elderly_dg1, fn_out8, row.names = FALSE)



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
