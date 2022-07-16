###########################################################
# spatial_DG_lifespan project
# Plot DEG for GO term "Regulation of Neurogenesis" results
# Anthony Ramnauth, June 01 2022
###########################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(pheatmap)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(sessioninfo)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on age bin
spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Infant")]

# Set gene names as row names for easier plotting
rownames(spe_pseudo) <- rowData(spe_pseudo)$gene_name

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data","pseudobulk_spe","NOinfant_modeling_results.rds"))

# Vector of GO "Regulation of Neurogenesis" term genes
GO_genes = as.character(c("PLXNA4", "SEMA3C", "SEMA5A", "TIAM1", "ULK1", "EPHA7",
    "SEMA5B", "CX3CL1", "DBN1", "PARP6", "EFNB3", "BAIAP2", "NPTN", "BDNF", "NUMBL",
    "NDEL1", "MBD1", "XRCC5", "DAB1", "CDK5R1", "SS18L2", "CAMK2B", "PRKCI", "ZFYVE27",
    "HDAC2", "SEMA6C", "NTRK3", "PPP3CA", "SORL1", "NEURL1", "BHLHB9", "SS18L1", "GDI1",
    "IL34", "ISLR2", "ITPKA", "PLXNA1", "SEMA7A", "NIN", "RTN4R", "MAP2K1", "CXCL12",
    "SEMA4F", "STK25", "RNF112", "FBXO31", "RTN4", "RGS14", "RNF10", "PRPF19", "ADNP",
    "TMEM98", "DUSP10", "DYNLT1", "DAG1", "TREM2", "MT3", "REST", "CXCR4", "EGFR",
    "CERS2", "PLXNB1", "MAP3K13", "CTNNA1", "LIMK1", "BHLHE40", "TNFRSF1B", "ASCL1",
    "FN1", "RELN", "PLXND1", "SIRT2", "DLX1", "SEMA3B", "NF2", "ASPA", "HES6",
    "NOTCH1", "GPER1", "SERPINF1", "NTRK2", "VEGFA", "SPP1", "LYN", "GJC2", "ETV5",
    "SEMA6A", "NKX6-2", "PTPRD", "TGM2", "HMGB2", "LRP4", "PAX6", "TSPO", "DAAM2",
    "IL6ST", "MAG", "CTDSP1", "OLIG2", "METRN", "SERPINE2", "SOX10", "PTN", "YAP1",
    "GPR37L1", "SOX8", "B2M", "KIT", "RND2", "BMP7"))

# Add logcounts for all clusters from GO term genes
neurogen_heatmap <- assays(spe_pseudo)[[2]][GO_genes,]
colnames(neurogen_heatmap) = paste("logcount", 1:48, sep = "")

# Add annotations for pheatmap
cluster_labels <- as.vector(c(rep("Cluster_1", 6), rep("Cluster_2", 6), rep("GCL", 6), rep("SGZ", 6),
    rep("CA4", 6), rep("CA3", 6), rep("ML", 6), rep("Cluster_8", 6)))

annotation_col <- data.frame(BayesSpace = factor(c(cluster_labels)))
rownames(annotation_col) = colnames(neurogen_heatmap)
ann_colors = list(BayesSpace = brewer.pal(8, "Set1"))
names(ann_colors$BayesSpace) <- unique(annotation_col$BayesSpace)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots","pseudobulked","NOinfant_neurogenesis_regulation_enrichment_heatmap_all.pdf"), width = 14, height = 14)
pheatmap(
    neurogen_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    color = inferno(20),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    fontsize_row = 9,
    main = "logcounts from enrichment model excluding infants with genes for neurogenesis regulation per cluster",
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

