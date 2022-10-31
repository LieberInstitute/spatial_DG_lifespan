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

# Load parent GO lists
load(file = here::here("processed-data", "pseudobulk_spe", "gene_ontologies", "DG_comp_parentGOdfs.Rdata"))

# Keep running into errors when trying to iteration a function (think because of 1:many mappings or invalid)
# so individually checking interesting GO terms

# Cellular Component Genes
inner_mito_membrane_protein_complex <- bitr("GO:0098800", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
respirasome <- bitr("GO:0070469", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
resp_chain_complex <- bitr("GO:0098803", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db) # invalid?
large_ribo_sub <- bitr("GO:0015934", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
lytic_vacuole_membrane <- bitr("GO:0098852", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db) # invaild?
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
intrinsic_comp_of_postsyn_membrane <- bitr("GO:0098936", fromType="GO",
    toType = "SYMBOL", OrgDb=org.Hs.eg.db) # invaild?
acetyltransferase_complex <- bitr("GO:1902493", fromType="GO",
    toType = "SYMBOL", OrgDb=org.Hs.eg.db) # invalid?

#######################################################
# Plot the GO parent genes for pseudobulk spe logcounts
#######################################################

respirasome <- unique(respirasome$SYMBOL)
respirasome <- respirasome[respirasome %in% rownames(spe_pseudo)]

# Add logcounts for all clusters from top age genes
respirasome <- assays(spe_pseudo)[[2]][respirasome, ]

Heatmap(respirasome,
    name = "logcounts",
    top_annotation = HeatmapAnnotation(age = spe_pseudo$age_bin, cluster = spe_pseudo$BayesSpace,
    col = list(age = c("Infant" = "purple", "Teen" = "blue", "Adult" = "red", "Elderly" = "forestgreen"),
        cluster = c("1" = "red4", "2" = "orange", "4" = "cyan", "8" = "springgreen3"))),
    column_title = "logcounts from GO:respirasome for GCL",
    show_column_names = FALSE,
    show_row_names = TRUE
    )



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
