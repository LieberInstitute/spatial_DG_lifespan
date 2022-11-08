#################################################################
# spatial_DG_lifespan project
# Plotting genes from parent CC GO terms of GCR across age groups
# Anthony Ramnauth, Oct 31 2022
#################################################################

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
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("1", "4", "8")]

# Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Load GO compare cluster
load(file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "GCR_comp_enrichedGO.Rdata"))

# Load parent GO lists
load(file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "GCR_comp_parentGOdfs.Rdata"))

#################################################
# Plot the GO parent terms in comparison dot plot
#################################################

CC_list <- CC_parent_GO$parent_term
MF_list <- MF_parent_GO$parent_term
BP_list <- BP_parent_GO$parent_term

pdf(file = here::here("plots", "pseudobulked", "Age_group_GCR_parent_GO.pdf"), width = 14, height = 12)

dotplot(comp_CC, showCategory = CC_list, label_format = 30, font.size = 20, includeAll = TRUE) +
    ggtitle("Parent GO terms for Cellular Compartment") +
    theme(plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = -45))

dotplot(comp_MF, showCategory = MF_list, label_format = 30, font.size = 20, includeAll = TRUE) +
    ggtitle("Parent GO terms for Molecular Function") +
    theme(plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = -45))

dotplot(comp_BP, showCategory = BP_list, label_format = 90, font.size = 14, includeAll = TRUE) +
    ggtitle("Parent GO terms for Biological Process")+
    theme(plot.title = element_text(size = 20),
        axis.text.x = element_text(angle = -45))

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC,
    showCategory = CC_list, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Parent GO terms for Cellular Compartment for Age groups in Granular Cell Region")

comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF,
    showCategory = MF_list, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Parent GO terms for Molecular Function for Age groups in Granular Cell Region")

comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP,
    showCategory = BP_list, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Parent GO terms for Biological Process for Age groups in Granular Cell Region")

dev.off()

#######################################################
# Plot the GO parent genes for pseudobulk spe logcounts
#######################################################

# Make function for heatmaps of all GO terms

heat<- function(x, y){

Heatmap(x,
    name = "logcounts",
    top_annotation = HeatmapAnnotation(age = spe_pseudo$age_bin, cluster = spe_pseudo$BayesSpace,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen"),
        cluster = c("1" = "red4", "4" = "cyan", "8" = "springgreen3"))),
    column_title = y,
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

CC_GOs <- CC_GOs[-15]

CC_GOs <- lapply(CC_GOs, findGOs)

# Add names to list

CC_names <- CC_parent_GO$parent_term
CC_names <- CC_names[-15]

names(CC_GOs) <- CC_names

CC_GOs <- lapply(CC_GOs, function(x) x%>% select(SYMBOL))

# Remove duplicates

CC_gene_list <- sapply(CC_GOs, unique)

# Find genes only in my data

CC_gene_list <- lapply(CC_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

# Plot

pdf(file = here::here("plots", "pseudobulked", "Parent_CC_GO_genes_GCR.pdf"), width = 14, height = 14)

mapply(heat, CC_gene_list, names(CC_gene_list))

dev.off()

# Molecular Function Genes

MF_GOs <- MF_parent_GO$parent_GO

# findGOs did not work for last GO term so remove and re-run

MF_GOs <- MF_GOs[-11]
MF_GOs <- MF_GOs[-11]

MF_GOs <- lapply(MF_GOs, findGOs)

# Add names to list

MF_parent_names <- MF_parent_GO$parent_term
MF_parent_names <- MF_parent_names[-11]
MF_parent_names <- MF_parent_names[-11]
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

pdf(file = here::here("plots", "pseudobulked", "Parent_MF_GO_genes_GCR.pdf"), width = 14, height = 14)

mapply(heat, MF_gene_list, names(MF_gene_list))

dev.off()

# Biological Processes Genes

BP_GOs <- BP_parent_GO$parent_GO

# findGOs did not work for last GO term so remove and re-run

BP_GOs <- BP_GOs[-10]
BP_GOs <- BP_GOs[-10]
BP_GOs <- BP_GOs[-16]
BP_GOs <- BP_GOs[-17]
BP_GOs <- BP_GOs[-17]
BP_GOs <- BP_GOs[-20]
BP_GOs <- BP_GOs[-23]
BP_GOs <- BP_GOs[-26]
BP_GOs <- BP_GOs[-27]
BP_GOs <- BP_GOs[-35]
BP_GOs <- BP_GOs[-39]

BP_GOs <- lapply(BP_GOs, findGOs)

# Add names to list

BP_names <- BP_parent_GO$parent_term
BP_names <- BP_names[-10]
BP_names <- BP_names[-10]
BP_names <- BP_names[-16]
BP_names <- BP_names[-17]
BP_names <- BP_names[-17]
BP_names <- BP_names[-20]
BP_names <- BP_names[-23]
BP_names <- BP_names[-26]
BP_names <- BP_names[-27]
BP_names <- BP_names[-35]
BP_names <- BP_names[-39]

names(BP_GOs) <- BP_names

BP_GOs <- lapply(BP_GOs, function(x) x%>% select(SYMBOL))

# Remove duplicates
BP_gene_list <- sapply(BP_GOs, unique)

# Find genes only in my data

BP_gene_list <- lapply(BP_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

# Remove GO terms that have errors in spe matrices (found by str(BP_gene_list))

BP_gene_list <- BP_gene_list[-16]
BP_gene_list <- BP_gene_list[-28]

pdf(file = here::here("plots", "pseudobulked", "Parent_BP_GO_genes_GCR.pdf"), width = 14, height = 14)

mapply(heat, BP_gene_list, names(BP_gene_list))

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
