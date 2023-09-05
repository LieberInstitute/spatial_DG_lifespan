##############################################
# spatial_DG_lifespan project
# DE analysis of age_bin vs. infant
# Anthony Ramnauth, Sept 2 2023
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

# Add colData() for label entire tissue minus CP
spe_pseudo$ALL <- 1

# Format spe object for DE models
colData(spe_pseudo) <- colData(spe_pseudo)[, sort(c(
    "age",
    "age_bin",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells",
    "ALL"
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

################################################################################################################
# Model and run DE for each spe_pseudo object of infant vs age_bin
################################################################################################################

# Create spe for Teen and Infant groups
teen_spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Adult")]
teen_spe_pseudo <- teen_spe_pseudo[, !teen_spe_pseudo$age_bin %in% c("Elderly")]

teen_spe_pseudo$enrichment_teen <- 0
teen_spe_pseudo$enrichment_teen[teen_spe_pseudo$age_bin == "Teen"] <- 1

model_formula <- ~enrichment_teen
m <- model.matrix(model_formula, data = colData(teen_spe_pseudo))

# Use pseudoBulkDGE to quickly perform DE analysis for entire tissue minus CP

teen_de_results <- pseudoBulkDGE(
    teen_spe_pseudo,
    col.data = colData(teen_spe_pseudo),
    label = teen_spe_pseudo$ALL,
    design = ~enrichment_teen,
    coef = "enrichment_teen",
    row.data = rowData(teen_spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
    )

# Create spe for Adult and Infant groups
adult_spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Teen")]
adult_spe_pseudo <- adult_spe_pseudo[, !adult_spe_pseudo$age_bin %in% c("Elderly")]

adult_spe_pseudo$enrichment_adult <- 0
adult_spe_pseudo$enrichment_adult[adult_spe_pseudo$age_bin == "Adult"] <- 1

model_formula <- ~enrichment_adult
m <- model.matrix(model_formula, data = colData(adult_spe_pseudo))

# Use pseudoBulkDGE to quickly perform DE analysis for entire tissue minus CP

adult_de_results <- pseudoBulkDGE(
    adult_spe_pseudo,
    col.data = colData(adult_spe_pseudo),
    label = adult_spe_pseudo$ALL,
    design = ~enrichment_adult,
    coef = "enrichment_adult",
    row.data = rowData(adult_spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
    )

# Create spe for Elderly and Infant groups
elderly_spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Teen")]
elderly_spe_pseudo <- elderly_spe_pseudo[, !elderly_spe_pseudo$age_bin %in% c("Adult")]

elderly_spe_pseudo$enrichment_elderly <- 0
elderly_spe_pseudo$enrichment_elderly[elderly_spe_pseudo$age_bin == "Elderly"] <- 1

model_formula <- ~enrichment_elderly
m <- model.matrix(model_formula, data = colData(elderly_spe_pseudo))

# Use pseudoBulkDGE to quickly perform DE analysis for entire tissue minus CP

elderly_de_results <- pseudoBulkDGE(
    elderly_spe_pseudo,
    col.data = colData(elderly_spe_pseudo),
    label = elderly_spe_pseudo$ALL,
    design = ~enrichment_elderly,
    coef = "enrichment_elderly",
    row.data = rowData(elderly_spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
    )

# Save modeling results
save(teen_de_results, adult_de_results, elderly_de_results,
    file = here::here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_age_bin_vs_infant_results.Rdata")
)

##########################
# Volcano Plots of results
##########################

DE_teen <- data.frame(
    gene_name = teen_de_results[[1]]$gene_name,
    logFC = teen_de_results[[1]]$logFC,
    adj.P.Val = teen_de_results[[1]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_teen <- ifelse(
    (DE_teen$adj.P.Val < 0.05) &
        (DE_teen$logFC > 1.5 | DE_teen$logFC < -1.5), "red", "gray47"
)

DE_teen_italics <- paste0("italic('", DE_teen$gene_name, "')")

## Legend names
names(keyvals_teen)[keyvals_teen == "red"] <- "Adjusted P-value < 0.05 & LFC > 1.5"
names(keyvals_teen)[keyvals_teen == "gray47"] <- "Not significant"

DE_adult <- data.frame(
    gene_name = adult_de_results[[1]]$gene_name,
    logFC = adult_de_results[[1]]$logFC,
    adj.P.Val = adult_de_results[[1]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_adult <- ifelse(
    (DE_adult$adj.P.Val < 0.05) &
        (DE_adult$logFC > 1.5 | DE_adult$logFC < -1.5), "red", "gray47"
)

DE_adult_italics <- paste0("italic('", DE_adult$gene_name, "')")

## Legend names
names(keyvals_adult)[keyvals_adult == "red"] <- "Adjusted P-value < 0.05 & LFC > 1.5"
names(keyvals_adult)[keyvals_adult == "gray47"] <- "Not significant"

DE_elderly <- data.frame(
    gene_name = elderly_de_results[[1]]$gene_name,
    logFC = elderly_de_results[[1]]$logFC,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_elderly <- ifelse(
    (DE_elderly$adj.P.Val < 0.05) &
        (DE_elderly$logFC > 1.5 | DE_elderly$logFC < -1.5), "red", "gray47"
)

DE_elderly_italics <- paste0("italic('", DE_elderly$gene_name, "')")

## Legend names
names(keyvals_elderly)[keyvals_elderly == "red"] <- "Adjusted P-value < 0.05 & LFC > 1.5"
names(keyvals_elderly)[keyvals_elderly == "gray47"] <- "Not significant"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_plots_age_bin_vs_infant.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(DE_teen,
    lab = DE_teen_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1.5,
    pCutoff = 0.049,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    colCustom = keyvals_teen,
    ylab = "-log10 Adjusted P-value",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "HPC",
    subtitle = "Teen vs. Infant",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-7, 7)) +
    ylim(c(0, 30))

EnhancedVolcano(DE_adult,
    lab = DE_adult_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1.5,
    pCutoff = 0.049,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    colCustom = keyvals_adult,
    ylab = "-log10 Adjusted P-value",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "HPC",
    subtitle = "Adult vs. Infant",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-7, 7)) +
    ylim(c(0, 30))

EnhancedVolcano(DE_elderly,
    lab = DE_elderly_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1.5,
    pCutoff = 0.049,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    colCustom = keyvals_elderly,
    ylab = "-log10 Adjusted P-value",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "HPC",
    subtitle = "Elderly vs. Infant",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-7, 7)) +
    ylim(c(0, 30))

dev.off()

######################################
# Write csv files for each DE analysis
######################################

# directory to save whole tissue results
dir_outputs <- here("processed-data", "pseudobulk_spe", "pseudoBulkDGE_results")

teen_infant <- data.frame(
    gene_id = teen_de_results[[1]]$gene_id,
    gene_name = teen_de_results[[1]]$gene_name,
    gene_type = teen_de_results[[1]]$gene_type,
    pvalue = teen_de_results[[1]]$P.Value,
    adj.P.Val = teen_de_results[[1]]$adj.P.Val,
    logFC = teen_de_results[[1]]$logFC
)

teen_infant <- teen_infant %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out1 <- file.path(dir_outputs, "Teen_vs_Infant_DE")

# Export summary as .csv file
write.csv(teen_infant, fn_out1, row.names = FALSE)

adult_infant <- data.frame(
    gene_id = adult_de_results[[1]]$gene_id,
    gene_name = adult_de_results[[1]]$gene_name,
    gene_type = adult_de_results[[1]]$gene_type,
    pvalue = adult_de_results[[1]]$P.Value,
    adj.P.Val = adult_de_results[[1]]$adj.P.Val,
    logFC = adult_de_results[[1]]$logFC
)

adult_infant <- adult_infant %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out2 <- file.path(dir_outputs, "Adult_vs_Infant_DE")

# Export summary as .csv file
write.csv(adult_infant, fn_out2, row.names = FALSE)

elderly_infant <- data.frame(
    gene_id = elderly_de_results[[1]]$gene_id,
    gene_name = elderly_de_results[[1]]$gene_name,
    gene_type = elderly_de_results[[1]]$gene_type,
    pvalue = elderly_de_results[[1]]$P.Value,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val,
    logFC = elderly_de_results[[1]]$logFC
)

elderly_infant <- elderly_infant %>%
    filter(adj.P.Val < 0.05) %>%
    dplyr::arrange(pvalue)

fn_out3 <- file.path(dir_outputs, "Elderly_vs_Infant_DE")

# Export summary as .csv file
write.csv(elderly_infant, fn_out3, row.names = FALSE)

###############################################################################

# Compile age signature gene set from significant DE genes from at least 2 DE results

combo_1 <- teen_infant[teen_infant$gene_id %in% adult_infant$gene_id,]
combo_2 <- teen_infant[teen_infant$gene_id %in% elderly_infant$gene_id,]
combo_3 <- adult_infant[adult_infant$gene_id %in% elderly_infant$gene_id,]

combo_1_up <- combo_1[combo_1$logFC >= 1.5,]
combo_1_up$sign <- 1
combo_1_down <- combo_1[combo_1$logFC <= -1.5,]
combo_1_down$sign <- -1
combo_2_up <- combo_2[combo_2$logFC >= 1.5,]
combo_2_up$sign <- 1
combo_2_down <- combo_2[combo_2$logFC <= -1.5,]
combo_2_down$sign <- -1
combo_3_up <- combo_3[combo_3$logFC >= 1.5,]
combo_3_up$sign <- 1
combo_3_down <- combo_3[combo_3$logFC <= -1.5,]
combo_3_down$sign <- -1

# Combine and find unique genes
combo_up <- rbind(combo_1_up, combo_2_up, combo_3_up)
combo_up <- unique(combo_up)

combo_down <- rbind(combo_1_down, combo_2_down, combo_3_down)
combo_down <- unique(combo_down)
