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

###################################
# Create sce from fetal HPC dataset
###################################

# Load fetal dataset directly into sce object

sce_fetal <- read10xCounts(here::here("GSE119212_hippocampus_aggr_8"))

rownames(sce_fetal) <- rowData(sce_fetal)$Symbol
colnames(sce_fetal) <- colData(sce_fetal)$Barcode
rownames(colData(sce_fetal)) <- colData(sce_fetal)$Barcode

dim(colData(sce_fetal))

# Add column data

fetal_coldata <- DataFrame(
    cell_ID = colData(sce_fetal)$Barcode,
    Dataset = rep("Zhong_etal_2020",times=33128)
)

colData(sce_fetal) <- fetal_coldata

# Save the fetal sce object

saveRDS(sce_fetal, file = here::here("sce_objects", "sce_fetal.rds"))

###############################################
# Create sce object for ant-to-post HPC dataset
###############################################

# Load ant to post axis data into a sparse matrix

atopaxis_sparse_mat <- readSparseCounts(here::here("GSE160189", "GSE160189_Hippo_Counts.csv.gz"), sep = ",",
    row.names = TRUE, col.names = TRUE)

dim(atopaxis_sparse_mat)

# Make temporary sce object from ant-to-post HPC dataset

sce_atopaxis <- SingleCellExperiment(assays = list(counts = atopaxis_sparse_mat))

dim(sce_atopaxis)

# Add column data

atopaxis_coldata <- DataFrame(
    cell_ID = colnames(sce_atopaxis),
    Dataset = rep("Ayhan_etal_2021",times=131325)
)

colData(sce_atopaxis) <- atopaxis_coldata
rownames(colData(sce_atopaxis)) <- sce_atopaxis$cell_ID

# Save the a to p axis sce data

saveRDS(sce_atopaxis, file = here::here("sce_objects", "sce_atopaxis.rds"))

##############################################################
# Compile Song lab datasets and combine into one sparse matrix
##############################################################

# One or some of these libraries have errors in resulting rownames causing duplication.
# Checking each file individually before merging to remove "\" character in rownames

song1_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5609934_sample1_humanHippocampus_1_lib1.dge.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

song1_2_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5609934_sample1_humanHippocampus_1_lib2.dge.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

song_red <- RowMergeSparseMatrices(song1_sparse_mat, song1_2_sparse_mat)

song9_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5610939_GEO_Sample9_seq1BC79_S4.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

song_red <- RowMergeSparseMatrices(song_red, song9_sparse_mat)

song13_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5610940_GEO_sample13_HCT15HBW_seq12BC82_S7.deg.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

song_red <- RowMergeSparseMatrices(song_red, song13_sparse_mat)

song15_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5610941_sample15_BC80_S5.deg.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

song_red <- RowMergeSparseMatrices(song_red, song15_sparse_mat)

song20_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5610942_Sample20_1.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song20_sparse_mat)

rownames(song20_sparse_mat) <- gsub("\"","",as.character(rownames(song20_sparse_mat)))
colnames(song20_sparse_mat) <- gsub("\"","",as.character(colnames(song20_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song20_sparse_mat)

song20_2_sparse_mat <- readSparseCounts(here::here("GSE185277_RAW",
    "GSM5610942_sample20_2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song20_2_sparse_mat)

rownames(song20_2_sparse_mat) <- gsub("\"","",as.character(rownames(song20_2_sparse_mat)))
colnames(song20_2_sparse_mat) <- gsub("\"","",as.character(colnames(song20_2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song20_2_sparse_mat)

song7_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618237_GEO_sample7_4yrs_plate1_seq12BC80_S5.deg.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song7_sparse_mat)

song_red <- RowMergeSparseMatrices(song_red, song7_sparse_mat)

song7_2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618237_GEO_sample7_4yrs_plate2_seq12BC81_S6.deg.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song7_2_sparse_mat)

song_red <- RowMergeSparseMatrices(song_red, song7_2_sparse_mat)

song6_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample6_B1_lib1.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song6_sparse_mat)
colnames(song6_sparse_mat)

rownames(song6_sparse_mat) <- gsub("\"","",as.character(rownames(song6_sparse_mat)))
colnames(song6_sparse_mat) <- gsub("\"","",as.character(colnames(song6_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song6_sparse_mat)

song6_2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample6-_B1_lib2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song6_2_sparse_mat)
colnames(song6_2_sparse_mat)

rownames(song6_2_sparse_mat) <- gsub("\"","",as.character(rownames(song6_2_sparse_mat)))
colnames(song6_2_sparse_mat) <- gsub("\"","",as.character(colnames(song6_2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song6_2_sparse_mat)

song11_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample11_B1_lib1.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song11_sparse_mat)

rownames(song11_sparse_mat) <- gsub("\"","",as.character(rownames(song11_sparse_mat)))
colnames(song11_sparse_mat) <- gsub("\"","",as.character(colnames(song11_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song11_sparse_mat)

song11_2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample11_B1_lib2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song11_2_sparse_mat)

rownames(song11_2_sparse_mat) <- gsub("\"","",as.character(rownames(song11_2_sparse_mat)))
colnames(song11_2_sparse_mat) <- gsub("\"","",as.character(colnames(song11_2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song11_2_sparse_mat)

song16_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample16_B1_lib1.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song16_sparse_mat)

rownames(song16_sparse_mat) <- gsub("\"","",as.character(rownames(song16_sparse_mat)))
colnames(song16_sparse_mat) <- gsub("\"","",as.character(colnames(song16_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song16_sparse_mat)

song16_2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample16_B1_lib2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song16_2_sparse_mat)

rownames(song16_2_sparse_mat) <- gsub("\"","",as.character(rownames(song16_2_sparse_mat)))
colnames(song16_2_sparse_mat) <- gsub("\"","",as.character(colnames(song16_2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song16_2_sparse_mat)

song22_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample22_B1_lib1.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song22_sparse_mat)

rownames(song22_sparse_mat) <- gsub("\"","",as.character(rownames(song22_sparse_mat)))
colnames(song22_sparse_mat) <- gsub("\"","",as.character(colnames(song22_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song22_sparse_mat)

song22_2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample22_B1_lib2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song22_2_sparse_mat)

rownames(song22_2_sparse_mat) <- gsub("\"","",as.character(rownames(song22_2_sparse_mat)))
colnames(song22_2_sparse_mat) <- gsub("\"","",as.character(colnames(song22_2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song22_2_sparse_mat)

song39_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample39_B1_lib1.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song39_sparse_mat)

rownames(song39_sparse_mat) <- gsub("\"","",as.character(rownames(song39_sparse_mat)))
colnames(song39_sparse_mat) <- gsub("\"","",as.character(colnames(song39_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song39_sparse_mat)

song39_2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618238_Sample39-lib2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song39_2_sparse_mat) <- gsub("\"","",as.character(rownames(song39_2_sparse_mat)))
colnames(song39_2_sparse_mat) <- gsub("\"","",as.character(colnames(song39_2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song39_2_sparse_mat)

song2_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618239_Sample2.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song2_sparse_mat) <- gsub("\"","",as.character(rownames(song2_sparse_mat)))
colnames(song2_sparse_mat) <- gsub("\"","",as.character(colnames(song2_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song2_sparse_mat)

song5_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618239_Sample5.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song5_sparse_mat) <- gsub("\"","",as.character(rownames(song5_sparse_mat)))
colnames(song5_sparse_mat) <- gsub("\"","",as.character(colnames(song5_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song5_sparse_mat)

song10_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618239_Sample10.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song10_sparse_mat) <- gsub("\"","",as.character(rownames(song10_sparse_mat)))
colnames(song10_sparse_mat) <- gsub("\"","",as.character(colnames(song10_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song10_sparse_mat)

song14_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618239_Sample14.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song14_sparse_mat) <- gsub("\"","",as.character(rownames(song14_sparse_mat)))
colnames(song14_sparse_mat) <- gsub("\"","",as.character(colnames(song14_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song14_sparse_mat)

song21_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618239_Sample21.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song21_sparse_mat) <- gsub("\"","",as.character(rownames(song21_sparse_mat)))
colnames(song21_sparse_mat) <- gsub("\"","",as.character(colnames(song21_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song21_sparse_mat)

song12_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618240_Sample12.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song12_sparse_mat) <- gsub("\"","",as.character(rownames(song12_sparse_mat)))
colnames(song12_sparse_mat) <- gsub("\"","",as.character(colnames(song12_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song12_sparse_mat)

song17_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618240_Sample17.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song17_sparse_mat) <- gsub("\"","",as.character(rownames(song17_sparse_mat)))
colnames(song17_sparse_mat) <- gsub("\"","",as.character(colnames(song17_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song17_sparse_mat)

song18_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618240_Sample18.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song18_sparse_mat) <- gsub("\"","",as.character(rownames(song18_sparse_mat)))
colnames(song18_sparse_mat) <- gsub("\"","",as.character(colnames(song18_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song18_sparse_mat)

song40_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618240_Sample40.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song40_sparse_mat) <- gsub("\"","",as.character(rownames(song40_sparse_mat)))
colnames(song40_sparse_mat) <- gsub("\"","",as.character(colnames(song40_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song40_sparse_mat)

song41_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618240_Sample41.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song41_sparse_mat) <- gsub("\"","",as.character(rownames(song41_sparse_mat)))
colnames(song41_sparse_mat) <- gsub("\"","",as.character(colnames(song41_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song41_sparse_mat)

song3_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618241_Sample3.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song3_sparse_mat) <- gsub("\"","",as.character(rownames(song3_sparse_mat)))
colnames(song3_sparse_mat) <- gsub("\"","",as.character(colnames(song3_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song3_sparse_mat)

song4_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618241_Sample4.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song4_sparse_mat) <- gsub("\"","",as.character(rownames(song4_sparse_mat)))
colnames(song4_sparse_mat) <- gsub("\"","",as.character(colnames(song4_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song4_sparse_mat)

song8_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618241_Sample8.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song8_sparse_mat) <- gsub("\"","",as.character(rownames(song8_sparse_mat)))
colnames(song8_sparse_mat) <- gsub("\"","",as.character(colnames(song8_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song8_sparse_mat)

song19_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618241_Sample19.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song19_sparse_mat) <- gsub("\"","",as.character(rownames(song19_sparse_mat)))
colnames(song19_sparse_mat) <- gsub("\"","",as.character(colnames(song19_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song19_sparse_mat)

song42_sparse_mat <- readSparseCounts(here::here("GSE185553_RAW",
    "GSM5618241_Sample42.txt.gz"), sep = "", row.names = TRUE, col.names = TRUE)

rownames(song42_sparse_mat) <- gsub("\"","",as.character(rownames(song42_sparse_mat)))
colnames(song42_sparse_mat) <- gsub("\"","",as.character(colnames(song42_sparse_mat)))

song_red <- RowMergeSparseMatrices(song_red, song42_sparse_mat)

# Load merged sparse matrix into an sce object

sce_song <- SingleCellExperiment(assays = list(counts = song_red))


# Add column data

dim(sce_song)

song_coldata <- DataFrame(
    cell_ID = colnames(sce_song),
    Dataset = rep("Zhou_etal_2022",times=187255)
)

stopifnot(identical(song_coldata$cell_ID, colnames(song_red)))

colnames(sce_song) <- sce_song$cell_ID

colData(sce_song) <- song_coldata

# Save the song lab sce data

saveRDS(sce_song, file = here::here("sce_objects", "sce_song.rds"))

####################
# Merge sce datasets
####################

# Not merging the a to p axis dataset because the number of features is very low

merged_sce_assays <- RowMergeSparseMatrices(assays(sce_sestan_DG)$counts, assays(sce_song)$counts)
merged_sce_assays <- RowMergeSparseMatrices(merged_sce_assays, assays(sce_fetal)$counts)

sce_sum <- SingleCellExperiment(assays = list(counts = merged_sce_assays))

# Merge coldata, but make sure you have cell_ID and  Dataset, in that order before rbind!

colData(sce_sestan_DG) <- colData(sce_sestan_DG)[, sort(c(
     "cell_ID",
     "Dataset"
 ))]

colData(sce_song) <- colData(sce_song)[, sort(c(
     "cell_ID",
     "Dataset"
 ))]

colData(sce_fetal) <- colData(sce_fetal)[, sort(c(
     "cell_ID",
     "Dataset"
 ))]

sum_coldata <- rbind(colData(sce_sestan_DG), colData(sce_song), colData(sce_fetal))

stopifnot(identical(rownames(sum_coldata), colnames(sce_sum)))

colData(sce_sum) <- sum_coldata

rowData(sce_sum)$SYMBOL <- rownames(sce_sum)

# Save the summed sce data

saveRDS(sce_sum, file = here::here("sce_objects", "sce_sum.rds"))
