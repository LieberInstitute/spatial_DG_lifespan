################################
# spatial_DG_lifespan project
# Sestan lab SCE obj creation
# Anthony Ramnauth, May 01 2023
################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scuttle)
    library(DropletUtils)
    library(here)
    library(scater)
    library(uwot)
    library(rtracklayer)
    library(Seurat)
})

#######################################
# First working with Sestan lab dataset
#######################################

# Load the sparse matrix

sestan_sparse_mat <- Matrix::readMM(here::here("GSE186538", "GSE186538_Human_counts.mtx.gz"))

dim(sestan_sparse_mat)

# Load gene list and metadata

sestan_genes <- read.delim(here::here("GSE186538", "GSE186538_Human_genes.txt.gz"), header = FALSE)

dim(sestan_genes)

sestan_ensembl <- read.delim(here::here("GSE186538", "ENSEMBL_ID_genes.tsv"), header = TRUE)

dim(sestan_ensembl)

sestan_meta_info <- read.delim(here::here("GSE186538", "GSE186538_Human_cell_meta.txt.gz"), sep = "\t", header = TRUE)

dim(sestan_meta_info)

# Add the rownames & colnames to the sparse matrix

rownames(sestan_sparse_mat) <- sestan_ensembl$ENSEMBL_ID
colnames(sestan_sparse_mat) <- sestan_meta_info$cell_name

# Make temporary sce object from sestan lab data

sce_sestan <- SingleCellExperiment(assays = list(counts = sestan_sparse_mat))

# Add coldata with metadata from sestan lab data

dim(sce_sestan)

sestan_coldata <- DataFrame(
    cell_name = sestan_meta_info[, "cell_name"],
    sample_ID = sestan_meta_info[, "samplename"],
    sum_umi = sestan_meta_info[, "nCount_RNA"],
    sum_gene = sestan_meta_info[, "nFeature_RNA"],
    mito_percent = sestan_meta_info[, "percent.mt"],
    region = sestan_meta_info[, "region"],
    cluster = sestan_meta_info[, "cluster"]
)

sestan_rowdata <- DataFrame(
    gene_name = sestan_ensembl$FINAL_NAME
)

stopifnot(identical(sestan_coldata$cell_name, colnames(sestan_sparse_mat)))

colData(sce_sestan) <- sestan_coldata
rowData(sce_sestan) <- sestan_rowdata
rownames(colData(sce_sestan)) <- sestan_coldata$cell_name

sce_sestan_DG <- sce_sestan[, sce_sestan$region == "DG"]

dim(colData(sce_sestan_DG))
dim(sce_sestan_DG)

# Remove any residual entorhinal cortex nuclei
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L5 BCL11B ADRA1A")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L6 TLE4 SULF1")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L6b TLE4 CCN2")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L3 PCP4 ADARB2")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L2 CUX2 CALB1")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L2 CUX2 IL1RAPL2")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L2 RELN BCL11B")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L5 RORB TLL1")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("EC L5 RORB TLL1")]

dim(sce_sestan_DG)
unique(sce_sestan_DG$cluster)

# Remove any residual subiculum nuclei
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("SUB distal FN1 NTNG1")]
sce_sestan_DG <- sce_sestan_DG[, !sce_sestan_DG$cluster %in% c("SUB proximal ROBO1 COL5A2")]

dim(sce_sestan_DG)
unique(sce_sestan_DG$cluster)

# Save the DG sestan lab sce object

saveRDS(sce_sestan_DG, file = here::here("processed-data", "sce", "sce_sestan_DG.rds"))
