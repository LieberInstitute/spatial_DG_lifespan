###############################
# spatial_DG_lifespan project
# Practice SCE obj creation
# Anthony Ramnauth, Nov 04 2022
###############################

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

sestan_meta_info <- read.delim(here::here("GSE186538", "GSE186538_Human_cell_meta.txt.gz"), sep = "\t", header = TRUE)

dim(sestan_meta_info)

# Add the rownames & colnames to the sparse matrix

rownames(sestan_sparse_mat) <- sestan_genes$V1
colnames(sestan_sparse_mat) <- sestan_meta_info$cell_name

# Make temporary sce object from sestan lab data

sce_sestan <- SingleCellExperiment(assays = list(counts = sestan_sparse_mat))

# Add coldata with metadata from sestan lab data

dim(sce_sestan)

sestan_coldata <- DataFrame(
    cell_ID = sestan_meta_info[, "cell_name"],
    Dataset = rep("Franjic_etal_2022",times=219058)
)

stopifnot(identical(sestan_coldata$cell_ID, colnames(sestan_sparse_mat)))

colData(sce_sestan) <- sestan_coldata
rownames(colData(sce_sestan)) <- sestan_coldata$cell_ID

sce_sestan_DG <- sce_sestan[, sce_sestan$region == "DG"]

dim(colData(sce_sestan_DG))

# Save the DG sestan lab sce object

saveRDS(sce_sestan_DG, file = here::here("sce_objects", "sce_sestan_DG.rds"))
