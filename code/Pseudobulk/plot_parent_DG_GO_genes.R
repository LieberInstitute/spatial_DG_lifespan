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
spe_pseudo <- spe_pseudo[, spe_pseudo$BayesSpace %in% c("1", "2", "4", "8")]

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

pdf(file = here::here("plots", "pseudobulked", "Age_group_DG_parent_GO.pdf"), width = 28, height = 21)

dotplot(comp_CC, showCategory = CC_list, label_format = 90, font.size = 26) +
    ggtitle("Parent GO terms for Cellular Compartment for Age groups in Dentate Gyrus") +
    theme(plot.title = element_text(size = 26),
        axis.text.x = element_text(angle = -45))

dotplot(comp_MF, showCategory = MF_list, label_format = 90, font.size = 26) +
    ggtitle("Parent GO terms for Molecular Function for Age groups in Dentate Gyrus") +
    theme(plot.title = element_text(size = 26),
        axis.text.x = element_text(angle = -45))

dotplot(comp_BP, showCategory = BP_list, label_format = 90, font.size = 26) +
    ggtitle("Parent GO terms for Biological Process for Age groups in Dentate Gyrus")+
    theme(plot.title = element_text(size = 26),
        axis.text.x = element_text(angle = -45))

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC,
    showCategory = CC_list, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Parent GO terms for Cellular Compartment for Age groups in Dentate Gyrus")

comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF,
    showCategory = MF_list, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Parent GO terms for Molecular Function for Age groups in Dentate Gyrus")

comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP,
    showCategory = BP_list, color = "p.adjust",
    pie = "count", cex_category = 3, label_format = 20
) +
    ggtitle("Parent GO terms for Biological Process for Age groups in Dentate Gyrus")

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
        cluster = c("1" = "red4", "2" = "orange", "4" = "cyan", "8" = "springgreen3"))),
    column_title = y,
    show_column_names = FALSE,
    show_row_names = TRUE
    )

}

# Keep running into errors when trying to iteration bitr function (think because of 1:many mappings or invalid)
# so individually checking interesting GO terms

# Cellular Component Genes
inner_mito_membrane_protein_complex <- bitr("GO:0098800", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
respirasome <- bitr("GO:0070469", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
large_ribo_sub <- bitr("GO:0015934", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
transport_vesicle <- bitr("GO:0030133", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
polysome <- bitr("GO:0005844", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
r_endo_ret_membrane <- bitr("GO:0030867", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
lysosomal_lumen <- bitr("GO:0043202", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
tertiary_granule <- bitr("GO:0070820", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
proteasome_complex <- bitr("GO:0000502", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
myelin_sheath <- bitr("GO:0043209", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
main_axon <- bitr("GO:0044304", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
ER_to_Golgi_transport_vesicle_membrane <- bitr("GO:0012507", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
cytochrome_complex <- bitr("GO:0070069", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
lumenal_side_of_membrane <- bitr("GO:0098576", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
postsynaptic_density <- bitr("GO:0014069", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
chaperone_complex <- bitr("GO:0101031", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
glutamatergic_synapse <- bitr("GO:0098978", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
RNA_poly_II_transcr_reg_complex <- bitr("GO:0090575", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
focal_adhesion <- bitr("GO:0005925", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
microtubule <- bitr("GO:0005874", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
cation_channel_complex <- bitr("GO:0034703", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)

CC_gene_list <- list(
    respirasome = unique(respirasome$SYMBOL),
    large_ribo_sub = unique(large_ribo_sub$SYMBOL),
    transport_vesicle = unique(transport_vesicle$SYMBOL),
    polysome = unique(polysome$SYMBOL),
    r_endo_ret_membrane = unique(r_endo_ret_membrane$SYMBOL),
    lysosomal_lumen = unique(lysosomal_lumen$SYMBOL),
    tertiary_granule = unique(tertiary_granule$SYMBOL),
    proteasome_complex = unique(proteasome_complex$SYMBOL),
    myelin_sheath = unique(myelin_sheath$SYMBOL),
    main_axon = unique(main_axon$SYMBOL),
    ER_to_Golgi_transport_vesicle_membrane = unique(ER_to_Golgi_transport_vesicle_membrane$SYMBOL),
    cytochrome_complex = unique(cytochrome_complex$SYMBOL),
    postsynaptic_density = unique(postsynaptic_density$SYMBOL),
    chaperone_complex = unique(chaperone_complex$SYMBOL),
    glutamatergic_synapse= unique(glutamatergic_synapse$SYMBOL),
    RNA_poly_II_transcr_reg_complex = unique(RNA_poly_II_transcr_reg_complex$SYMBOL),
    focal_adhesion= unique(focal_adhesion$SYMBOL),
    microtubule = unique(microtubule$SYMBOL),
    cation_channel_complex = unique(cation_channel_complex$SYMBOL)
)


CC_gene_list <- lapply(CC_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

pdf(file = here::here("plots", "pseudobulked", "Parent_CC_GO_genes_DentateGyrus.pdf"), width = 14, height = 14)

mapply(heat, CC_gene_list, names(CC_gene_list))

dev.off()

# Molecular Function Genes

MF_GOs <- MF_parent_GO$parent_GO

# findGOs did not work for last GO term so remove

MF_GOs <- MF_GOs[-15]

findGOs <- function(x){
    bitr(x, fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
}

MF_GOs <- lapply(MF_GOs, findGOs)

# Add names to list

MF_parent_names <- MF_parent_GO$parent_term[-15]
names(MF_GOs) <- MF_parent_names

MF_GOs <- lapply(MF_GOs, function(x) x%>% select(SYMBOL))

MF_gene_list <- sapply(MF_GOs, unique)

# Remove duplicate


MF_gene_list <- lapply(MF_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

pdf(file = here::here("plots", "pseudobulked", "Parent_MF_GO_genes_DentateGyrus.pdf"), width = 14, height = 14)

mapply(heat, MF_gene_list, names(MF_gene_list))

dev.off()

# Biological Processes Genes

BP_GOs <- BP_parent_GO$parent_GO

# findGOs did not work for last GO term so remove

BP_GOs <- BP_GOs[-15]
BP_GOs <- BP_GOs[-20]
BP_GOs <- BP_GOs[-32]
BP_GOs <- BP_GOs[-36]
BP_GOs <- BP_GOs[-36]
BP_GOs <- BP_GOs[-41]
BP_GOs <- BP_GOs[-41]
BP_GOs <- BP_GOs[-46]
BP_GOs <- BP_GOs[-46]
BP_GOs <- BP_GOs[-46]
BP_GOs <- BP_GOs[-49]
BP_GOs <- BP_GOs[-51]
BP_GOs <- BP_GOs[-53]

BP_GOs <- lapply(BP_GOs, findGOs)

# Add names to list

BP_names <- BP_parent_GO$parent_term
BP_names <- BP_names[-15]
BP_names <- BP_names[-20]
BP_names <- BP_names[-32]
BP_names <- BP_names[-36]
BP_names <- BP_names[-36]
BP_names <- BP_names[-41]
BP_names <- BP_names[-41]
BP_names <- BP_names[-46]
BP_names <- BP_names[-46]
BP_names <- BP_names[-46]
BP_names <- BP_names[-49]
BP_names <- BP_names[-51]
BP_names <- BP_names[-53]


BP_parent_names <- BP_names
names(BP_GOs) <- BP_parent_names

BP_GOs <- lapply(BP_GOs, function(x) x%>% select(SYMBOL))

BP_gene_list <- sapply(BP_GOs, unique)

# Remove duplicate


BP_gene_list <- lapply(BP_gene_list, function(x) {

    ix <- x[x %in% rownames(spe_pseudo)]

    assays(spe_pseudo)[[2]][ix, ]

})

BP_gene_list <- BP_gene_list[-41]
BP_gene_list <- BP_gene_list[-30]
BP_gene_list <- BP_gene_list[-39]
BP_gene_list <- BP_gene_list[-39]
BP_gene_list <- BP_gene_list[-41]
BP_gene_list <- BP_gene_list[-42]

pdf(file = here::here("plots", "pseudobulked", "Parent_BP_GO_genes_DentateGyrus.pdf"), width = 14, height = 14)

mapply(heat, BP_gene_list, names(BP_gene_list))

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
