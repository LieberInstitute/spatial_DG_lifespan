###############################
# spatial_DG_lifespan project
# Hochgerner SCE obj creation
# Anthony Ramnauth, Jan 26 2024
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
    library(Matrix)
})

# load the raw data

hochgerner <- read.delim(here::here("GSE104323_10X_expression_data_V2.tab.gz"), sep = "\t", header = TRUE)

# Need to make some adjustments to get a sparse matrix
rownames(hochgerner) <- hochgerner$cellid

hochgerner$cellid <- NULL

hochgerner1 <- as.matrix(hochgerner)

hochgerner2 <- as(hochgerner1, "dgCMatrix")

# Load metadata
hochgerner_meta <- read.delim(here::here("GSE104323_metadata_barcodes_24185cells.txt.gz"), sep = "\t", header = TRUE)

# Make temporary sce object from hochgerner data

sce <- SingleCellExperiment(assays = list(counts = hochgerner2))

# Add metadata to colData(sce)

hochgerner_coldata <- DataFrame(
    cell_name = hochgerner_meta[, "Sample.name..24185.single.cells."],
    strain = hochgerner_meta[, "characteristics..strain"],
    age = hochgerner_meta[, "characteristics..age"],
    pooled = hochgerner_meta[, "characteristics..sex.of.pooled.animals"],
    cellType = hochgerner_meta[, "characteristics..cell.cluster"],
    molecule = hochgerner_meta[, "molecule"],
    SRR.run.accession = hochgerner_meta[, "SRR.run.accession"],
    raw.file.name = hochgerner_meta[, "raw.file..original.file.name."],
    UMI_CellularBarcode = hochgerner_meta[, "UMI_CellularBarcode"]
)

rownames(hochgerner_coldata) <- hochgerner_coldata$cell_name

stopifnot(identical(rownames(hochgerner_coldata), colnames(sce)))
# Error: identical(hochgerner_coldata$cell_name, colnames(hochgerner2)) is not TRUE
# Noticed a difference in colnames and rownames - and .
hochgerner_coldata$cell_name2 <- sub("-$", ".", hochgerner_coldata$cell_name)
rownames(hochgerner_coldata) <- hochgerner_coldata$cell_name2
rownames(hochgerner_coldata) <- paste0("X", rownames(hochgerner_coldata))

# remove rows whose rownames don't match the assay(sce, "counts")
matching_rows <- rownames(hochgerner_coldata) %in% colnames(sce)

# Subset the dataframe based on the matching rows
hochgerner_coldata_subset <- hochgerner_coldata[matching_rows, ]

stopifnot(identical(rownames(hochgerner_coldata_subset), colnames(sce)))

# match the ordering
match_rows <- match(colnames(sce), rownames(hochgerner_coldata_subset))

# Reorder the rows in the dataframe
hochgerner_coldata_subset <- hochgerner_coldata_subset[match_rows, ]

stopifnot(identical(rownames(hochgerner_coldata_subset), colnames(sce)))

colData(sce) <- hochgerner_coldata_subset

rowData(sce)$gene_name <- rownames(sce)

# Save the Hochgerner HPC  sce object

saveRDS(sce, file = here::here("sce_objects", "sce_hochgerner.rds"))

