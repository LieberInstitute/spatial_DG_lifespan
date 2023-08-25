###########################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked dentate gyrus across age_bin
# Anthony Ramnauth, Sept 26 2022
###########################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

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
spe_pseudo$dentate_gyrus[spe_pseudo$BayesSpace == "2" |
        spe_pseudo$BayesSpace == "4"|
        spe_pseudo$BayesSpace == "6"|
        spe_pseudo$BayesSpace == "7"] <- 1

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

# Use pseudoBulkDGE to quickly perform age_bin DE analysis for BayesSpace labels
infant_de_results <- pseudoBulkDGE(
    spe_pseudo,
    col.data = colData(spe_pseudo),
    label = spe_pseudo$dentate_gyrus,
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
    label = spe_pseudo$dentate_gyrus,
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
    label = spe_pseudo$dentate_gyrus,
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
    label = spe_pseudo$dentate_gyrus,
    design = ~enrichment_elderly,
    coef = "enrichment_elderly",
    row.data = rowData(spe_pseudo),
    method = "voom",
    qualities = TRUE,
    robust = TRUE
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
    adj.P.Val = infant_de_results[[1]]$adj.P.Val
)

dentate1infant <- data.frame(
    gene_name = infant_de_results[[2]]$gene_name,
    logFC = infant_de_results[[2]]$logFC,
    adj.P.Val = infant_de_results[[2]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_inf <- ifelse(
    dentate1infant$adj.P.Val < 0.05, "darksalmon", "#f0e3d6"
)

# Create groups of genes according to their function

Neurogenic <- paste0(
    "italic('",
    c("SOX11",
	"NES",
	"DCX",
	"DPYSL5",
	"BHLHE22",
    "DIO2"),
    "')")

Activated_microglia <- paste0(
    "italic('",
    c(
    "CD74",
	"C1QC",
	"C1QA",
	"HLA-DQB1",
	"HLA-DQA1",
    "HLA−DPB1",
    "HLA−DPA1",
    "HLA−DRA",
    "TREM2"),
    "')")

Reactive_astro <- paste0(
    "italic('",
    c(
    "GFAP",
	"CD44",
	"C3",
	"CD14",
	"SERPINA3"),
    "')")

GABAergic <- paste0(
    "italic('",
    c(
    "GAD1",
    "GAD2",
    "LAMP5",
    "CCK",
    "VIP",
    "SST",
    "PVALB",
    "RELN",
    "TAC1",
    "NPY"
    ),
    "')")

selected <- c(Neurogenic, Activated_microglia, Reactive_astro, GABAergic)

dentate1infant_italics <- paste0("italic('", dentate1infant$gene_name, "')")

## Assigning colors for each groups of highlited genes
keyvals_inf[dentate1infant_italics %in% Activated_microglia] <- "#A5C0DF"
keyvals_inf[dentate1infant_italics %in% Neurogenic] <- "#789C25"
keyvals_inf[dentate1infant_italics %in% Reactive_astro] <- "#006164"
keyvals_inf[dentate1infant_italics %in% GABAergic] <- "orange"

## Legend names
names(keyvals_inf)[keyvals_inf == "darksalmon"] <- "Adjusted P-value < 0.05"
names(keyvals_inf)[keyvals_inf == "#f0e3d6"] <- "Not significant"
names(keyvals_inf)[keyvals_inf == "#789C25"] <- "Neurogenic genes"
names(keyvals_inf)[keyvals_inf == "#A5C0DF"] <- "Activated microglia genes"
names(keyvals_inf)[keyvals_inf == "#006164"] <- "Reactive astroglia genes"
names(keyvals_inf)[keyvals_inf == "orange"] <- "GABAergic genes"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Infant_dentategyrus.pdf"),
    width = 10.5, height = 8)

EnhancedVolcano(dentate1infant,
    lab = dentate1infant_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = selected,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    pointSize = c(ifelse(dentate1infant_italics %in% selected, 6, 2)),
    colAlpha = c(ifelse(dentate1infant_italics %in% selected, 1, 0.2)),
    colCustom = keyvals_inf,
    max.overlaps = Inf,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Infant vs. non-Infant",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-3.5, 3.5)) +
    ylim(c(0, 14))

dev.off()

###########################################
# Volcano Plots of results for teen age_bin
###########################################

dentate0teen <- data.frame(
    gene_name = teen_de_results[[1]]$gene_name,
    logFC = teen_de_results[[1]]$logFC,
    adj.P.Val = teen_de_results[[1]]$adj.P.Val
)

dentate1teen <- data.frame(
    gene_name = teen_de_results[[2]]$gene_name,
    logFC = teen_de_results[[2]]$logFC,
    adj.P.Val = teen_de_results[[2]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_teen <- ifelse(
    dentate1teen$adj.P.Val < 0.05, "darksalmon", "#f0e3d6"
)

dentate1teen_italics <- paste0("italic('", dentate1teen$gene_name, "')")

## Assigning colors for each groups of highlited genes
keyvals_teen[dentate1teen_italics %in% Activated_microglia] <- "#A5C0DF"
keyvals_teen[dentate1teen_italics %in% Neurogenic] <- "#789C25"
keyvals_teen[dentate1teen_italics %in% Reactive_astro] <- "#006164"
keyvals_teen[dentate1teen_italics %in% GABAergic] <- "orange"

## Legend names
names(keyvals_teen)[keyvals_teen == "darksalmon"] <- "Adjusted P-value < 0.05"
names(keyvals_teen)[keyvals_teen == "#f0e3d6"] <- "Not significant"
names(keyvals_teen)[keyvals_teen == "#789C25"] <- "Neurogenic genes"
names(keyvals_teen)[keyvals_teen == "#A5C0DF"] <- "Activated microglia genes"
names(keyvals_teen)[keyvals_teen == "#006164"] <- "Reactive astroglia genes"
names(keyvals_teen)[keyvals_teen == "orange"] <- "GABAergic genes"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Teen_dentategyrus.pdf"),
    width = 10.5, height = 8)

EnhancedVolcano(dentate1teen,
    lab = dentate1teen_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = selected,
    parseLabels = TRUE,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    pointSize = c(ifelse(dentate1teen_italics %in% selected, 6, 2)),
    colAlpha = c(ifelse(dentate1teen_italics %in% selected, 1, 0.2)),
    colCustom = keyvals_teen,
    max.overlaps = Inf,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Teen vs. non-Teen",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-2.5, 1)) +
    ylim(c(0, 2.5))

dev.off()

############################################
# Volcano Plots of results for adult age_bin
############################################

dentate0adult <- data.frame(
    gene_name = adult_de_results[[1]]$gene_name,
    logFC = adult_de_results[[1]]$logFC,
    adj.P.Val = adult_de_results[[1]]$adj.P.Val
)

dentate1adult <- data.frame(
    gene_name = adult_de_results[[2]]$gene_name,
    logFC = adult_de_results[[2]]$logFC,
    adj.P.Val = adult_de_results[[2]]$adj.P.Val
)


## Colors for the significant and not significant genes
keyvals_adult <- ifelse(
    dentate1adult$adj.P.Val < 0.05, "darksalmon", "#f0e3d6"
)

dentate1adult_italics <- paste0("italic('", dentate1adult$gene_name, "')")

## Assigning colors for each groups of highlited genes
keyvals_adult[dentate1adult_italics %in% Activated_microglia] <- "#A5C0DF"
keyvals_adult[dentate1adult_italics %in% Neurogenic] <- "#789C25"
keyvals_adult[dentate1adult_italics %in% Reactive_astro] <- "#006164"
keyvals_adult[dentate1adult_italics %in% GABAergic] <- "orange"

## Legend names
names(keyvals_adult)[keyvals_adult == "darksalmon"] <- "Adjusted P-value < 0.05"
names(keyvals_adult)[keyvals_adult == "#f0e3d6"] <- "Not significant"
names(keyvals_adult)[keyvals_adult == "#789C25"] <- "Neurogenic genes"
names(keyvals_adult)[keyvals_adult == "#A5C0DF"] <- "Activated microglia genes"
names(keyvals_adult)[keyvals_adult == "#006164"] <- "Reactive astroglia genes"
names(keyvals_adult)[keyvals_adult == "orange"] <- "GABAergic genes"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Adult_dentategyrus.pdf"),
    width = 10.5, height = 8)

EnhancedVolcano(dentate1adult,
    lab = dentate1adult_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = selected,
    parseLabels = TRUE,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    pointSize = c(ifelse(dentate1adult_italics %in% selected, 6, 2)),
    colAlpha = c(ifelse(dentate1adult_italics %in% selected, 1, 0.2)),
    colCustom = keyvals_adult,
    max.overlaps = Inf,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Adult vs. non-Adult",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-1.5, 1.5)) +
    ylim(c(0, 3))

dev.off()

##############################################
# Volcano Plots of results for elderly age_bin
##############################################

dentate0elderly <- data.frame(
    gene_name = elderly_de_results[[1]]$gene_name,
    logFC = elderly_de_results[[1]]$logFC,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val
)

dentate1elderly <- data.frame(
    gene_name = elderly_de_results[[2]]$gene_name,
    logFC = elderly_de_results[[2]]$logFC,
    adj.P.Val = elderly_de_results[[2]]$adj.P.Val
)

## Colors for the significant and not significant genes
keyvals_elderly <- ifelse(
    dentate1elderly$adj.P.Val < 0.05, "darksalmon", "#f0e3d6"
)

dentate1elderly_italics <- paste0("italic('", dentate1elderly$gene_name, "')")

## Assigning colors for each groups of highlited genes
keyvals_elderly[dentate1elderly_italics %in% Activated_microglia] <- "#A5C0DF"
keyvals_elderly[dentate1elderly_italics %in% Neurogenic] <- "#789C25"
keyvals_elderly[dentate1elderly_italics %in% Reactive_astro] <- "#006164"
keyvals_elderly[dentate1elderly_italics %in% GABAergic] <- "orange"

## Legend names
names(keyvals_elderly)[keyvals_elderly == "darksalmon"] <- "Adjusted P-value < 0.05"
names(keyvals_elderly)[keyvals_elderly == "#f0e3d6"] <- "Not significant"
names(keyvals_elderly)[keyvals_elderly == "#789C25"] <- "Neurogenic genes"
names(keyvals_elderly)[keyvals_elderly == "#A5C0DF"] <- "Activated microglia genes"
names(keyvals_elderly)[keyvals_elderly == "#006164"] <- "Reactive astroglia genes"
names(keyvals_elderly)[keyvals_elderly == "orange"] <- "GABAergic genes"

pdf(file = here::here("plots", "pseudobulked","pseudoBulkDGE", "pseudoBulkDGE_DE_volcano_Elderly_dentategyrus.pdf"),
    width = 10.5, height = 8)

EnhancedVolcano(dentate1elderly,
    lab = dentate1elderly_italics,
    x = 'logFC',
    y = 'adj.P.Val',
    FCcutoff = 1,
    pCutoff = 0.049,
    selectLab = selected,
    parseLabels = TRUE,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    pointSize = c(ifelse(dentate1elderly_italics %in% selected, 6, 2)),
    colAlpha = c(ifelse(dentate1elderly_italics %in% selected, 1, 0.2)),
    colCustom = keyvals_elderly,
    max.overlaps = Inf,
    ylab = "-log10 adj.P.Val",
    legendLabels = c('Not sig.','Log (base 2) FC','adj.P.Val',
      'adj.P.Val & Log (base 2) FC'),
    title = "Dentate Gyrus",
    subtitle = "Elderly vs. non-Elderly",
    legendPosition = "bottom",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-2, 3.5)) +
    ylim(c(0, 13))

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
    pvalue = infant_de_results[[1]]$P.Value,
    adj.P.Val = infant_de_results[[1]]$adj.P.Val,
    logFC = infant_de_results[[1]]$logFC
)

infant_dg0 <- infant_dg0 %>%
    dplyr::arrange(pvalue)

fn_out1 <- file.path(dir_outputs, "InfantvsNonInfant_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(infant_dg0, fn_out1, row.names = FALSE)

infant_dg1 <- data.frame(
    gene_id = infant_de_results[[2]]$gene_id,
    gene_name = infant_de_results[[2]]$gene_name,
    gene_type = infant_de_results[[2]]$gene_type,
    pvalue = infant_de_results[[2]]$P.Value,
    adj.P.Val = infant_de_results[[2]]$adj.P.Val,
    logFC = infant_de_results[[2]]$logFC
)

infant_dg1 <- infant_dg1 %>%
    dplyr::arrange(pvalue)

fn_out2 <- file.path(dir_outputs, "InfantvsNonInfant_DentateGyrus_DE")

# Export summary as .csv file
write.csv(infant_dg1, fn_out2, row.names = FALSE)

teen_dg0 <- data.frame(
    gene_id = teen_de_results[[1]]$gene_id,
    gene_name = teen_de_results[[1]]$gene_name,
    gene_type = teen_de_results[[1]]$gene_type,
    pvalue = teen_de_results[[1]]$P.Value,
    adj.P.Val = teen_de_results[[1]]$adj.P.Val,
    logFC = teen_de_results[[1]]$logFC
)

teen_dg0 <- teen_dg0 %>%
    dplyr::arrange(pvalue)

fn_out3 <- file.path(dir_outputs, "TeenvsNonTeen_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(teen_dg0, fn_out3, row.names = FALSE)

teen_dg1 <- data.frame(
    gene_id = teen_de_results[[2]]$gene_id,
    gene_name = teen_de_results[[2]]$gene_name,
    gene_type = teen_de_results[[2]]$gene_type,
    pvalue = teen_de_results[[2]]$P.Value,
    adj.P.Val = teen_de_results[[2]]$adj.P.Val,
    logFC = teen_de_results[[2]]$logFC
)

teen_dg1 <- teen_dg1 %>%
    dplyr::arrange(pvalue)

fn_out4 <- file.path(dir_outputs, "TeenvsNonTeen_DentateGyrus_DE")

# Export summary as .csv file
write.csv(teen_dg1, fn_out4, row.names = FALSE)

adult_dg0 <- data.frame(
    gene_id = adult_de_results[[1]]$gene_id,
    gene_name = adult_de_results[[1]]$gene_name,
    gene_type = adult_de_results[[1]]$gene_type,
    pvalue = adult_de_results[[1]]$P.Value,
    adj.P.Val = adult_de_results[[1]]$adj.P.Val,
    logFC = adult_de_results[[1]]$logFC
)

adult_dg0 <- adult_dg0 %>%
    dplyr::arrange(pvalue)

fn_out5 <- file.path(dir_outputs, "AdultvsNonAdult_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(adult_dg0, fn_out5, row.names = FALSE)

adult_dg1 <- data.frame(
    gene_id = adult_de_results[[2]]$gene_id,
    gene_name = adult_de_results[[2]]$gene_name,
    gene_type = adult_de_results[[2]]$gene_type,
    pvalue = adult_de_results[[2]]$P.Value,
    adj.P.Val = adult_de_results[[2]]$adj.P.Val,
    logFC = adult_de_results[[2]]$logFC
)

adult_dg1 <- adult_dg1 %>%
    dplyr::arrange(pvalue)

fn_out6 <- file.path(dir_outputs, "AdultvsNonAdult_DentateGyrus_DE")

# Export summary as .csv file
write.csv(adult_dg1, fn_out6, row.names = FALSE)

elderly_dg0 <- data.frame(
    gene_id = elderly_de_results[[1]]$gene_id,
    gene_name = elderly_de_results[[1]]$gene_name,
    gene_type = elderly_de_results[[1]]$gene_type,
    pvalue = elderly_de_results[[1]]$P.Value,
    adj.P.Val = elderly_de_results[[1]]$adj.P.Val,
    logFC = elderly_de_results[[1]]$logFC
)

elderly_dg0 <- elderly_dg0 %>%
    dplyr::arrange(pvalue)

fn_out7 <- file.path(dir_outputs, "ElderlyvsNonElderly_NonDentateGyrus_DE")

# Export summary as .csv file
write.csv(elderly_dg0, fn_out7, row.names = FALSE)

elderly_dg1 <- data.frame(
    gene_id = elderly_de_results[[2]]$gene_id,
    gene_name = elderly_de_results[[2]]$gene_name,
    gene_type = elderly_de_results[[2]]$gene_type,
    pvalue = elderly_de_results[[2]]$P.Value,
    adj.P.Val = elderly_de_results[[2]]$adj.P.Val,
    logFC = elderly_de_results[[2]]$logFC
)

elderly_dg1 <- elderly_dg1 %>%
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
