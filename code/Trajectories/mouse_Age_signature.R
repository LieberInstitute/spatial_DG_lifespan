#####################################################
# spatial_DG_lifespan project
# Mouse common Age signature analysis on human Visium
# Anthony Ramnauth, Oct 19 2023
#####################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(sessioninfo)
    library(SingleCellExperiment)
    library(rafalib)
    library(limma)
    library(edgeR)
    library(scran)
    library(VISION)
    library(dplyr)
    library(ggh4x)
    library(ggsignif)
    library(spatialLIBD)
    library(lsmeans)
    library(RColorBrewer)
})

# Load mouse CAS genes
CAS <- read.csv(file = here("processed-data","gene_set_enrichment",
    "Hahn_2023.csv"))

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Translate from one species to the other using the orthology
CAS <- orthology[orthology$Column3 %in% CAS$Gene.Symbol,]
CAS <- unique(CAS$Column1)

# Convert to ensembIDs that are present in my Visium dataset
# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

CAS <- rowData(spe)$gene_id[match(CAS, rowData(spe)$gene_name)]

# remove NA
CAS <- CAS[!is.na(CAS)]

# Make into dataframe for signed geneset
CAS <- data.frame(
    gene_id = CAS,
    sign = rep(1, 82)
)

values_to_change <- c("ENSG00000144230", "ENSG00000120694", "ENSG00000132002", "ENSG00000110172",
    "ENSG00000122884", "ENSG00000100591")

CAS <- CAS %>%
  mutate(sign = case_when(
    gene_id %in% values_to_change ~ -1,
    TRUE ~ sign
  ))

# Create a named vector or VISION gene signature creation

AgingSignature <- setNames(CAS$sign, CAS$gene_id)

########################
# VISION analysis of CAS
########################

sig <- createGeneSignature(name = "CommonAgingSignature", sigData = AgingSignature)

mySignatures <- c(sig)

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
  colData(spe)$mouse_CAS <- as.numeric(plyr::mapvalues(x = row.names(colData(spe)),
      from = names(sigScores), to = as.numeric(sigScores)))

# order spe observations according to age for better plotting
spe <- spe[, order(spe$age)]

# Plot CAS on tissue ordered by age
vis_grid_gene(
    spe = spe,
    geneid = "mouse_CAS",
    pdf = here::here("plots", "Trajectories", "mouse_CAS_BayesSpace_spotplot.pdf"),
    minCount = 0,
    viridis = FALSE,
    point_size = 1.5,
    spatial = TRUE,
    image_id = "lowres"
    )

scores <- spe$mouse_CAS

#load spe object that has human CAS for checking correlation of mosue and human scores on Visium spots
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_CAS_BayesSpace_spe.rds"))

# add to spe object
spe$mouse_CAS <- scores

# Save new spe object
saveRDS(spe,
    file = here::here("processed-data", "harmony_processed_spe", "harmony_CAS_BayesSpace_spe.rds"))

# Load SPE with cell type mean abundances
spe_cell <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

# order spe observations according to age
spe_cell <- spe_cell[, order(spe_cell$age)]

# Make data frame with  cell abundances and use key as common rownames
colnames(spe_cell) <- spe_cell$key

cell_abundances <- as.data.frame(colData(spe_cell)[, 44:68])

# Get dominant cell types per Visium spot for plots

cell_props <- as.data.frame(colData(spe_cell)[, c(44:68)])
for (col in 1:ncol(cell_props)){
    colnames(cell_props)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(cell_props)[col])
}

# Find dominant celltype for each spot
p <- cell_props %>%
 rowwise() %>%
 mutate(row_max = names(.)[which.max(c_across(everything()))])

colData(spe_cell)$dominant_cell_types <- p$row_max

# Set colors
cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3",
    "VLMC" = "firebrick1", "Pericyte" = "red2", "Endoth" = "red", "SMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "magenta", "Astro_2" = "yellow3")

new_data <- data.frame(
    human_CAS = spe$CAS,
    mouse_CAS = spe$mouse_CAS
)

# Fit a linear model to data
fit_cas <- lm(human_CAS ~ mouse_CAS, data = new_data)

# Extract the coefficients and R2 value
intercept <- coef(fit_cas)[1]
slope <- coef(fit_cas)[2]
r2 <- summary(fit_cas)$r.squared

pdf(file = here::here("plots", "Trajectories", "Human_CAS_vs_Mouse_CAS.pdf"), width = 9, height = 5)

ggplot(new_data,
    aes(
        x = mouse_CAS,
        y = human_CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("Mouse CAS") +
    ylab("Human CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = -0.3, y = 2.3,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept, r2),
        size = 4) +
    annotate("text", x = -0.5, y = 2.1,
        label = sprintf("R2 = %.2f", r2), size = 4)

dev.off()
