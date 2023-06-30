################################################################
# spatial_DG_lifespan project
# Plotting genes from parent CC GO terms of DG across age groups
# Anthony Ramnauth, Oct 31 2022
################################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on BayesSpace clusters for DG
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("2", "4", "6", "7")]

bayes_df <- data.frame(spe_pseudo$BayesSpace)
bayes_df <- bayes_df %>%
    mutate(DG_layer = case_when(
        grepl("2", spe_pseudo.BayesSpace) ~ "ML",
        grepl("4", spe_pseudo.BayesSpace) ~ "CA3&4",
        grepl("6", spe_pseudo.BayesSpace) ~ "SGZ",
        grepl("7", spe_pseudo.BayesSpace) ~ "GCL"
    ))

colData(spe_pseudo)$BayesSpace <- factor(bayes_df$DG_layer, levels = c("ML", "CA3&4", "SGZ", "GCL"))

# Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Load GO compare cluster
load(file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "DG_comp_enrichedGO.Rdata"))

# Load parent GO lists
load(file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "DG_comp_parentGOdfs.Rdata"))

#################################################
# Plot the GO parent terms in comparison dot plot
#################################################

CC_list <- CC_parent_GO$parent_term
MF_list <- MF_parent_GO$parent_term
BP_list <- BP_parent_GO$parent_term

# Truncate the BP_list, it is too large for plotting

trunc_CC_list <- c(
    "adherens junction",
    "histone methyltransferase complex",
    "lamellipodium",
    "apical junction complex",
    "NuRD complex",
    "postsynaptic density membrane",
    "methyltransferase complex",
    "lysosomal membrane",
    "transport vesicle",
    "mitochondrial matrix",
    "vacuolar proton−transporting V-type ATPase complex",
    "proteasome complex",
    "MHC class II protein complex",
    "endocytic vesicle",
    "cytoplasmic stress granule",
    "translation preinitiation complex",
    "respirasome",
    "presynapse",
    "stress fiber"
)

trunc_MF_list <- c(
    "transcription coregulator activity",
    "histone binding",
    "histone demethylase activity",
    "ATP-dependent chromatin remodeler activity",
    "miRNA binding",
    "lyase activity",
    "MHC class II protein complex binding",
    "ATPase−coupled ion transmembrane transporter activity",
    "electron transfer activity",
    "translation initiation factor activity",
    "oxidoreduction−driven active transmembrane transporter activity",
    "antioxidant activity",
    "tau-protein kinase activity",
    "immune receptor activity",
    "MHC class II receptor activity"
)

trunc_BP_list <- c(
    "regulation of neuron projection development",
    "regulation of nervous system development",
    "regulation of neurogenesis",
    "synapse organization",
    "histone modification",
    "gliogenesis",
    "neural precursor cell proliferation",
    "axon development",
    "dendrite development",
    "dendritic spine development",
    "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II",
    "aerobic respiration",
    "myeloid cell activation involved in immune response",
    "regulation of translation",
    "cellular respiration",
    "leukocyte mediated immunity",
    "fear response",
    "'de novo' post-translational protein folding",
    "maintenance of blood-brain barrier"
)

pdf(file = here::here("plots", "pseudobulked", "Select_Age_group_DG_parent_GO.pdf"), width = 8, height = 9)

dotplot(comp_CC, showCategory = trunc_CC_list, label_format = 90) +
    ggtitle("Parent GO terms for Cellular Compartment for Age groups in Dentate Gyrus") +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70, size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

dotplot(comp_MF, showCategory = trunc_MF_list, label_format = 90) +
    ggtitle("Parent GO terms for Molecular Function for Age groups in Dentate Gyrus") +
    theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = -70),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

dotplot(comp_BP, showCategory = trunc_BP_list) +
    ggtitle("Select Parent GO terms for Biological Process for Age groups in Dentate Gyrus")+
    theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(angle = -70),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

dev.off()

# Plot dotplots separately to format sizes

pdf(file = here::here("plots", "pseudobulked", "blah_Age_group_DG_parent_GOBP_separate.pdf"), width = 20, height = 32)

dotplot(comp_BP, showCategory = BP_list, label_format = 300) +
    ggtitle("GO terms for CC for Age groups in Dentate Gyrus")+
    theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(angle = -70, size = 16),
        axis.text.y = element_text(size = 16, face = "bold"))

dev.off()


#######################################################
# Plot the GO parent genes for pseudobulk spe logcounts
#######################################################

# Make function for heatmaps of all GO terms

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
    16, 12, 13, 15, 8, 1, 11, 2, 9, 10, 14, 4, 3, 5, 7, 6,
    32, 28, 29, 31, 24, 17, 27, 18, 25, 26, 30, 20, 19, 21, 23, 22,
    48, 44, 45, 47, 40, 33, 43, 34, 41, 42, 46, 36, 35, 37, 39, 38,
    64, 60, 61, 63, 56, 49, 59, 50, 57, 58, 62, 52, 51, 53, 55, 54
)

# Make function for heatmaps of all GO terms

heat<- function(x, y){

Heatmap(x,
    name = "logcounts",
    top_annotation = HeatmapAnnotation(cluster = spe_pseudo$BayesSpace, age = spe_pseudo$age_bin,
    col = list(cluster = c("ML" = "#E4E1E3", "CA3&4" = "#FEAF16", "SGZ" = "#1CFFCE", "GCL" = "#B00068"),
        age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen")
        )),
    column_title = y,
    column_order = Bayes_age_order,
    show_column_names = FALSE,
    show_row_names = TRUE
    )

}

# Keep running into errors when trying to iteration bitr function (think because of 1:many mappings or invalid)
# so individually checking interesting GO terms

# Cellular Component Genes

CC_GOs <- CC_parent_GO$parent_GO

findGOs <- function(x){
    bitr(x, fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
}

# findGOs did not work for last GO term so remove and re-run

CC_GOs <- lapply(CC_GOs, findGOs)

# Add names to list

CC_names <- CC_parent_GO$parent_term
names(CC_GOs) <- CC_names

CC_GOs <- lapply(CC_GOs, function(x) x%>% select(SYMBOL))

# Remove duplicates

CC_gene_list <- sapply(CC_GOs, unique)

# Find genes only in my data

CC_gene_list <- lapply(CC_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

# If needed, remove GO terms that have errors in spe matrices (found by str(CC_gene_list))
CC_gene_list <- CC_gene_list[-8]
CC_gene_list <- CC_gene_list[-11]

# Plot

pdf(file = here::here("plots", "pseudobulked", "Parent_CC_GO_genes_DentateGyrus.pdf"), width = 14, height = 14)

mapply(heat, CC_gene_list, names(CC_gene_list))

dev.off()

# Molecular Function Genes

MF_GOs <- MF_parent_GO$parent_GO

# findGOs did not work for last GO term so remove and re-run

MF_GOs <- MF_GOs[-2]
MF_GOs <- MF_GOs[-5]
MF_GOs <- MF_GOs[-19]

MF_GOs <- lapply(MF_GOs, findGOs)

# Add names to list (make sure to remove the same entries that failed earlier)

MF_parent_names <- MF_parent_GO$parent_term[-2]
MF_parent_names <- MF_parent_names[-5]
MF_parent_names <- MF_parent_names[-19]

names(MF_GOs) <- MF_parent_names

MF_GOs <- lapply(MF_GOs, function(x) x%>% select(SYMBOL))

# Remove duplicates

MF_gene_list <- sapply(MF_GOs, unique)

# Find genes only in my data

MF_gene_list <- lapply(MF_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

# Plot

pdf(file = here::here("plots", "pseudobulked", "Parent_MF_GO_genes_DentateGyrus.pdf"), width = 14, height = 14)

mapply(heat, MF_gene_list, names(MF_gene_list))

dev.off()

# Biological Processes Genes

BP_GOs <- BP_parent_GO$parent_GO

# findGOs did not work for last GO term so remove and re-run

BP_GOs <- BP_GOs[-3]
BP_GOs <- BP_GOs[-8]
BP_GOs <- BP_GOs[-30]
BP_GOs <- BP_GOs[-43]
BP_GOs <- BP_GOs[-64]
BP_GOs <- BP_GOs[-78]
BP_GOs <- BP_GOs[-79]
BP_GOs <- BP_GOs[-92]
BP_GOs <- BP_GOs[-92]
BP_GOs <- BP_GOs[-102]
BP_GOs <- BP_GOs[-102]
BP_GOs <- BP_GOs[-104]
BP_GOs <- BP_GOs[-105]
BP_GOs <- BP_GOs[-106]
BP_GOs <- BP_GOs[-106]
BP_GOs <- BP_GOs[-116]
BP_GOs <- BP_GOs[-122]
BP_GOs <- BP_GOs[-128]
BP_GOs <- BP_GOs[-128]
BP_GOs <- BP_GOs[-134]
BP_GOs <- BP_GOs[-141]
BP_GOs <- BP_GOs[-150]

BP_GOs <- lapply(BP_GOs, findGOs)

# Add names to list

BP_names <- BP_parent_GO$parent_term
BP_names <- BP_names[-3]
BP_names <- BP_names[-8]
BP_names <- BP_names[-30]
BP_names <- BP_names[-43]
BP_names <- BP_names[-64]
BP_names <- BP_names[-78]
BP_names <- BP_names[-79]
BP_names <- BP_names[-92]
BP_names <- BP_names[-92]
BP_names <- BP_names[-102]
BP_names <- BP_names[-102]
BP_names <- BP_names[-104]
BP_names <- BP_names[-105]
BP_names <- BP_names[-106]
BP_names <- BP_names[-106]
BP_names <- BP_names[-116]
BP_names <- BP_names[-122]
BP_names <- BP_names[-128]
BP_names <- BP_names[-128]
BP_names <- BP_names[-134]
BP_names <- BP_names[-141]
BP_names <- BP_names[-150]

BP_parent_names <- BP_names
names(BP_GOs) <- BP_parent_names

BP_GOs <- lapply(BP_GOs, function(x) x%>% select(SYMBOL))

# Remove duplicates
BP_gene_list <- sapply(BP_GOs, unique)

# Find genes only in my data

BP_gene_list <- lapply(BP_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

# Remove GO terms that have errors in spe matrices size (found by str(BP_gene_list))

BP_gene_list <- BP_gene_list[-4]
BP_gene_list <- BP_gene_list[-9]
BP_gene_list <- BP_gene_list[-46]
BP_gene_list <- BP_gene_list[-52]
BP_gene_list <- BP_gene_list[-52]
BP_gene_list <- BP_gene_list[-65]
BP_gene_list <- BP_gene_list[-70]
BP_gene_list <- BP_gene_list[-81]
BP_gene_list <- BP_gene_list[-81]
BP_gene_list <- BP_gene_list[-92]
BP_gene_list <- BP_gene_list[-97]
BP_gene_list <- BP_gene_list[-3]
BP_gene_list <- BP_gene_list[-23]
BP_gene_list <- BP_gene_list[-60]
BP_gene_list <- BP_gene_list[-110]
BP_gene_list <- BP_gene_list[-118]
BP_gene_list <- BP_gene_list[-118]
BP_gene_list <- BP_gene_list[-127]
BP_gene_list <- BP_gene_list[-131]


pdf(file = here::here("plots", "pseudobulked", "Parent_BP_GO_genes_DentateGyrus.pdf"), width = 14, height = 14)

mapply(heat, BP_gene_list, names(BP_gene_list))

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
