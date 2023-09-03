###############################
# spatial_DG_lifespan project
# Common Age signature analysis
# Anthony Ramnauth, Sept 2 2023
###############################

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
    library(VISION)
    library(dplyr)
    library(ggh4x)
    library(spatialLIBD)
    library(lsmeans)
    library(RColorBrewer)
})

# Load pseudo-bulked SPE
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

##################################################################################################
# Compile common age signature (CAS) gene set from significant DE genes from at least 2 DE results
##################################################################################################


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

# Make sure there are no intersect between up and down genes
intersect(combo_up$gene_name,combo_down$gene_name)

CAS <- rbind(combo_up, combo_down)
CAS$pvalue <- NULL
CAS$adj.P.Val <- NULL
CAS$logFC <- NULL
CAS <- unique(CAS)

fn_out4 <- file.path(dir_outputs, "CAS_list")

# Export summary as .csv file
write.csv(CAS, fn_out4, row.names = FALSE)

# Create a named vector or VISION gene signature creation

AgingSignature <- setNames(CAS$sign, CAS$gene_id)

########################
# VISION analysis of CAS
########################

sig <- createGeneSignature(name = "CommonAgingSignature", sigData = AgingSignature)

mySignatures <- c(sig)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# For unique labels use spe$key (some barcodes are duplicated as labels)
colnames(spe) <- spe$key

# Scale counts within a sample
n.umi <- colSums(assays(spe)$counts)

scaled_counts <- t(t(assays(spe)$counts) / n.umi) * median(n.umi)

vis <- Vision(scaled_counts,
              signatures = mySignatures)

# Set the number of threads when running parallel computations
# On Windows, this must either be omitted or set to 1
options(mc.cores = 1)

vis <- analyze(vis)

#The analysis has finished, and we'll remap now the scores for each of the signatures
  sigScores <- vis@SigScores[,1]
  sigRanks <- rank(sigScores)
  colData(spe)$CAS <- as.numeric(plyr::mapvalues(x = row.names(colData(spe)),
      from = names(sigScores), to = as.numeric(sigScores)))

# order spe observations according to age for better plotting
spe <- spe[, order(spe$age)]

# Plot CAS on tissue ordered by age
vis_grid_gene(
    spe = spe,
    geneid = "CAS",
    pdf = here::here("plots", "Trajectories", "CAS_spotplot.pdf"),
    minCount = 0,
    viridis = FALSE,
    alpha = 0.5,
    point_size = 2,
    spatial = TRUE,
    image_id = "lowres",
    auto_crop = TRUE
    )

# Setting up data for violin plots of CAS vs age_bin faceted by spatial domain

df <-
    data.frame(spe$key, spe$sample_id, spe$bayesSpace_harmony_10)

df <- df %>%
    mutate(
        BayesSpace = case_when(
            spe.bayesSpace_harmony_10 == 1 ~ "SLM",
            spe.bayesSpace_harmony_10 == 2 ~ "ML",
            spe.bayesSpace_harmony_10 == 3 ~ "CP",
            spe.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe.bayesSpace_harmony_10 == 5 ~ "SR",
            spe.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe.bayesSpace_harmony_10 == 7 ~ "GCL",
            spe.bayesSpace_harmony_10 == 8 ~ "SL",
            spe.bayesSpace_harmony_10 == 9 ~ "CA1",
            spe.bayesSpace_harmony_10 == 10 ~ "WM",
        )
    )

colData(spe)$BayesSpace <-
    factor(df$BayesSpace, levels = c("SLM", "ML", "CP", "CA3&4", "SR", "SGZ",
        "GCL", "SL", "CA1", "WM"))

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        between(spe.age, 0, 3) ~ "Infant",
        between(spe.age, 13, 19) ~ "Teen",
        between(spe.age, 20, 50) ~ "Adult",
        between(spe.age, 60, 100) ~ "Elderly"
    ))

stopifnot(age_df$spe.key == spe$key)

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

saveRDS(spe,
    file = here::here("processed-data", "harmony_processed_spe", "harmony_CAS_spe.rds"))

# Create datafame of cell proportions and DG layer
CAS_df <- as.data.frame(colData(spe)[, c(44:46)],
    row.names = spe$key)

# Add ages
CAS_df$age <- spe$age

# Remove CP cluster since that has higher variance and was not included in upstream DE analyses
# Remove Choroid Plexus cluster
CAS_df$spatialnum <- spe$bayesSpace_harmony_10
CAS_df = CAS_df[which(CAS_df$spatialnum != "3"), ]

strip <- strip_themed(background_x = elem_list_rect(fill = c("#5A5156", "#E4E1E3", "#FEAF16",
    "#FE00FA", "#1CFFCE", "#B00068", "#90AD1C", "#16FF32", "#2ED9FF")))

bay_colors <- c("#5A5156", "#E4E1E3", "#FEAF16",
    "#FE00FA", "#1CFFCE", "#B00068", "#90AD1C", "#16FF32", "#2ED9FF")

pdf(file = here::here("plots", "Trajectories", "CAS_vs_age.pdf"), width = 7, height = 5)

ggplot(CAS_df, aes(x = age, y = CAS)) +
    geom_boxplot(width = 10, notch = F, outlier.colour = NA,  mapping = aes(fill = as.factor(age))) +
    geom_smooth(data = CAS_df, mapping = aes(x = age, y = CAS), color='red3', se = F, method = 'loess') +
    facet_wrap2(~ BayesSpace, nrow = 2, ncol = 5, strip = strip) +
    theme(axis.text.x = element_text(size = 7), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank()) +
    guides(fill = FALSE) +
    ylim(-2, 1) +
    theme_classic()

ggplot(CAS_df, aes(x = CAS_df$age_bin, y = CAS_df$CAS)) +
    geom_violin(aes(fill = age_bin)) +
    scale_fill_manual(values=c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")) +
    geom_boxplot(width = 0.1) +
    geom_smooth(data = CAS_df, mapping = aes(x = as.numeric(age_bin), y = CAS), color='red3', se = F, method = 'loess') +
    labs(x = "age", y = "CAS") +
    ggtitle("Common Aging Signature") +
    theme_classic() +
    ylim(-2, 1) +
    facet_wrap2(~ BayesSpace, nrow = 2, ncol = 5, strip = strip) +
    theme(axis.text.x = element_text(size = 7), legend.position = "none", plot.title = element_blank(),
        axis.title.x = element_blank())

ggplot(CAS_df, aes(x = age, y = CAS, color = BayesSpace)) +
    geom_smooth(data = CAS_df, mapping = aes(x = age, y = CAS), se = F,  method = 'loess') +
    scale_color_manual(values = bay_colors) +
    scale_y_continuous(limits = c(-0.3, 0.5)) +
    theme_classic()

ggplot(CAS_df, aes(x = age, y = CAS, color = BayesSpace)) +
    geom_smooth(data = CAS_df, mapping = aes(x = age, y = CAS), se = F,  method = 'lm') +
    scale_color_manual(values = bay_colors) +
    scale_y_continuous(limits = c(-0.3, 0.5)) +
    theme_classic()

dev.off()

# Linear modeling for CAS velocities

# set the spatial domain information to a new, character-type column, as the factors can confuse the linear models
CAS_df$spatial_domain <- as.character(CAS_df$BayesSpace)

# Calculate a linear model that allows for an interaction between age and spatial domain
m_interaction <- lm(CAS ~ age*spatial_domain, data = CAS_df)

# Setup a dataframe that contains the intercepts for each tissue for later inspection of the offset
coefficient_dataframe <- as.data.frame(m_interaction$coefficients)
colnames(coefficient_dataframe) <- 'coefficients'
coefficient_dataframe <- coefficient_dataframe[grep('^spatial_domain', row.names(coefficient_dataframe), value = T),,drop=F]
coefficient_dataframe$spatial_domain <- gsub('spatial_domain', '',row.names(coefficient_dataframe))
colnames(coefficient_dataframe)[1] <- 'offset'

# Use the function lstrends (from the lsmeans package) to allow domain-to-domain comparison
m_lst <- lstrends(m_interaction, "spatial_domain", var="age")

#Create pairwise comparisons between slopes and test using the pairs function and store that result
#This lets us know which slopes are statistically significantly different from one another
slope_pair_wiseTestResults <- as.data.frame(pairs(m_lst))

# Store it in a separate dataframe to hold the information about the slopes
Visium_slope_dataframe <- as.data.frame(m_lst)

# Take the coefficient dataframe that we generate above together with the slope dataframe
Visium_coefficient_dataframe <- merge.data.frame(Visium_slope_dataframe, coefficient_dataframe, by = 'spatial_domain')

# Test correlation between Baseline CAS and slope
cor.test(Visium_coefficient_dataframe$offset, Visium_coefficient_dataframe$age.trend)

bay_colors_noCA1 <- c("SLM" = "#5A5156", "ML" = "#E4E1E3", "CA3&4" = "#FEAF16",
   "SR" = "#FE00FA", "SGZ" = "#1CFFCE", "GCL" = "#B00068", "SL" = "#90AD1C", "WM" = "#2ED9FF")

# Fit a linear model to the slope data for plot annotation
fit_trend <- lm(age.trend ~ offset, data = Visium_coefficient_dataframe)

# Extract the coefficients and R2 value
intercept_trend <- coef(fit_trend)[1]
slope_trend <- coef(fit_trend)[2]
r2_trend <- summary(fit_trend)$r.squared

pdf(file = here::here("plots", "Trajectories", "CAS_slope_vs_baseline.pdf"), width = 7, height = 5)

ggplot(Visium_coefficient_dataframe, aes(x = offset, y = age.trend, color = spatial_domain)) +
    scale_color_manual(values = bay_colors_noCA1) +
    geom_point()  +
    labs(x = "Offset (a.u.)", y = "CAS Slope (a.u.)") +
    geom_smooth(method = 'lm', color='black', linetype = 'dashed', se = F) +
    scale_y_continuous(limits = c(0,0.01))  +
    scale_x_continuous(limits = c(-0.1,0.15)) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black'),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 0.1, y = 0.0095,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_trend, r2_trend),
        size = 4) +
    annotate("text", x = 0.1, y = 0.0090,
        label = sprintf("R2 = %.2f", r2_trend), size = 4) +
    theme_classic()

dev.off()

# Compare CAS slope per spatial domain in barplot

pdf(file = here::here("plots", "Trajectories", "CAS_slope_spatial_domains.pdf"))

ggplot(Visium_coefficient_dataframe, aes(x = reorder(spatial_domain, -age.trend), y = age.trend, fill = age.trend)) +
    geom_bar(position = position_dodge(), stat = "identity", mapping = aes(fill = age.trend)) +
     xlab("") +
     ylab("CAS slope")  +
    labs(fill = "CAS slope") +
    scale_fill_gradient2(low = brewer.pal(11, "RdYlBu")[10], midpoint = quantile(Visium_coefficient_dataframe$age.trend)[3],
        mid = brewer.pal(11, "RdYlBu")[6], high = brewer.pal(11, "RdYlBu")[2])+
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = .2,                    # Width of the error bars
                  position  =position_dodge(.9)) +
     scale_y_continuous(expand = c(0,0)) +
     theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1)) +
    theme_classic()

dev.off()
