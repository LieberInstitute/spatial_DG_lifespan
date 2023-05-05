#########################################
# spatial_DG_lifespan project
# NeuronChat spot-spot communication
# Anthony Ramnauth, May 04 2023
#########################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(SpatialExperiment)
    library(NeuronChat)
    library(CellChat)
    library(patchwork)
    library(spatialLIBD)
})

options(stringsAsFactors = FALSE)

# Load SPE
spe <-
    readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Ensure that each barcode is unique use spe$key for colnames
colnames(spe) <- spe$key

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        between(spe.age, 0, 3) ~ "Infant",
        between(spe.age, 13, 19) ~ "Teen",
        between(spe.age, 20, 50) ~ "Adult",
        between(spe.age, 60, 100) ~ "Elderly"
    ))

stopifnot(age_df$spe.key == spe$key)

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

# Set rownames as gene names since ensemble IDs seem not to work
rownames(spe) <- rowData(spe)$gene_name

# For ease of referencing, use DG region names instead of BayesSpace cluster #
Bayes_df <-
    data.frame(spe$key, spe$sample_id, spe$bayesSpace_harmony_10)
Bayes_df <- Bayes_df %>%
    mutate(
        BayesSpace = case_when(
            spe.bayesSpace_harmony_10 == 2 ~ "ML",
            spe.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe.bayesSpace_harmony_10 == 7 ~ "GCL",
        )
    )

colData(spe)$BayesSpace <-
    factor(Bayes_df$BayesSpace, levels = c("ML", "CA3&4", "SGZ", "GCL"))

# Subset for age bin
spe_infant <- spe[, spe$age_bin == "Infant"]
spe_not_infant <- spe[, !spe$age_bin == "Infant"]
spe_teen <- spe[, spe$age_bin == "Teen"]
spe_not_teen <- spe[, !spe$age_bin == "Teen"]
spe_adult <- spe[, spe$age_bin == "Adult"]
spe_not_adult <- spe[, !spe$age_bin == "Adult"]
spe_elderly <- spe[, spe$age_bin == "Elderly"]
spe_not_elderly <- spe[, !spe$age_bin == "Elderly"]

# Subset for Dentate Gyrus BayesSpace clusters
###############################################################################
spe_infant <- spe_infant[, spe_infant$bayesSpace_harmony_10 == "2" |
        spe_infant$bayesSpace_harmony_10 == "4" |
        spe_infant$bayesSpace_harmony_10 == "6" |
        spe_infant$bayesSpace_harmony_10 == "7"]

infant_matrix <- as.matrix(assays(spe_infant)$logcounts)
stopifnot(colnames(infant_matrix) == rownames(colData(spe_infant)))

spe_not_infant <- spe_not_infant[, spe_not_infant$bayesSpace_harmony_10 == "2" |
        spe_not_infant$bayesSpace_harmony_10 == "4" |
        spe_not_infant$bayesSpace_harmony_10 == "6" |
        spe_not_infant$bayesSpace_harmony_10 == "7"]

non_infant_matrix <- as.matrix(assays(spe_not_infant)$logcounts)
stopifnot(colnames(non_infant_matrix) == rownames(colData(spe_not_infant)))
###############################################################################

spe_teen <- spe_teen[, spe_teen$bayesSpace_harmony_10 == "2" |
        spe_teen$bayesSpace_harmony_10 == "4" |
        spe_teen$bayesSpace_harmony_10 == "6" |
        spe_teen$bayesSpace_harmony_10 == "7"]

teen_matrix <- as.matrix(assays(spe_teen)$logcounts)
stopifnot(colnames(teen_matrix) == rownames(colData(spe_teen)))

spe_not_teen <- spe_not_teen[, spe_not_teen$bayesSpace_harmony_10 == "2" |
        spe_not_teen$bayesSpace_harmony_10 == "4" |
        spe_not_teen$bayesSpace_harmony_10 == "6" |
        spe_not_teen$bayesSpace_harmony_10 == "7"]

non_teen_matrix <- as.matrix(assays(spe_not_teen)$logcounts)
stopifnot(colnames(non_teen_matrix) == rownames(colData(spe_not_teen)))
###############################################################################

spe_adult <- spe_adult[, spe_adult$bayesSpace_harmony_10 == "2" |
        spe_adult$bayesSpace_harmony_10 == "4" |
        spe_adult$bayesSpace_harmony_10 == "6" |
        spe_adult$bayesSpace_harmony_10 == "7"]

adult_matrix <- as.matrix(assays(spe_adult)$logcounts)
stopifnot(colnames(adult_matrix) == rownames(colData(spe_adult)))

spe_not_adult <- spe_not_adult[, spe_not_adult$bayesSpace_harmony_10 == "2" |
        spe_not_adult$bayesSpace_harmony_10 == "4" |
        spe_not_adult$bayesSpace_harmony_10 == "6" |
        spe_not_adult$bayesSpace_harmony_10 == "7"]

non_adult_matrix <- as.matrix(assays(spe_not_adult)$logcounts)
stopifnot(colnames(non_adult_matrix) == rownames(colData(spe_not_adult)))
###############################################################################

spe_elderly <- spe_elderly[, spe_elderly$bayesSpace_harmony_10 == "2" |
            spe_elderly$bayesSpace_harmony_10 == "4" |
            spe_elderly$bayesSpace_harmony_10 == "6" |
            spe_elderly$bayesSpace_harmony_10 == "7"]

elderly_matrix <- as.matrix(assays(spe_elderly)$logcounts)
stopifnot(colnames(elderly_matrix) == rownames(colData(spe_elderly)))

spe_not_elderly <- spe_not_elderly[, spe_not_elderly$bayesSpace_harmony_10 == "2" |
            spe_not_elderly$bayesSpace_harmony_10 == "4" |
            spe_not_elderly$bayesSpace_harmony_10 == "6" |
            spe_not_elderly$bayesSpace_harmony_10 == "7"]

non_elderly_matrix <- as.matrix(assays(spe_not_elderly)$logcounts)
stopifnot(colnames(non_elderly_matrix) == rownames(colData(spe_not_elderly)))
###############################################################################

# Create neuronchat object for age binned spe objects
neuronchat_infant <- createNeuronChat(object = infant_matrix, DB='human',
    group.by = spe_infant$BayesSpace)

neuronchat_not_infant <- createNeuronChat(object = non_infant_matrix, DB='human',
    group.by = spe_not_infant$BayesSpace)

neuronchat_teen <- createNeuronChat(object = teen_matrix, DB='human',
    group.by = spe_teen$BayesSpace)

neuronchat_not_teen <- createNeuronChat(object = non_teen_matrix, DB='human',
    group.by = spe_not_teen$BayesSpace)

neuronchat_adult <- createNeuronChat(object = adult_matrix, DB='human',
    group.by = spe_adult$BayesSpace)

neuronchat_not_adult <- createNeuronChat(object = non_adult_matrix, DB='human',
    group.by = spe_not_adult$BayesSpace)

neuronchat_elderly <- createNeuronChat(object = elderly_matrix, DB='human',
    group.by = spe_elderly$BayesSpace)

neuronchat_not_elderly <- createNeuronChat(object = non_elderly_matrix, DB='human',
    group.by = spe_not_elderly$BayesSpace)

############################
# Run NeuronChat
############################

neuronchat_infant <- run_NeuronChat(neuronchat_infant,M=100)
neuronchat_not_infant <- run_NeuronChat(neuronchat_not_infant,M=100)
neuronchat_teen <- run_NeuronChat(neuronchat_teen,M=100)
neuronchat_not_teen <- run_NeuronChat(neuronchat_not_teen,M=100)
neuronchat_adult <- run_NeuronChat(neuronchat_adult,M=100)
neuronchat_not_adult <- run_NeuronChat(neuronchat_not_adult,M=100)
neuronchat_elderly <- run_NeuronChat(neuronchat_elderly,M=100)
neuronchat_not_elderly <- run_NeuronChat(neuronchat_not_elderly,M=100)

################
# Merge datasets
################

neuronchat_infant@meta <- as.data.frame(colData(spe_infant), row.names = rownames(colData(spe_infant)))
neuronchat_not_infant@meta <- as.data.frame(colData(spe_not_infant), row.names = rownames(colData(spe_not_infant)))
neuronchat_teen@meta <- as.data.frame(colData(spe_teen), row.names = rownames(colData(spe_teen)))
neuronchat_not_teen@meta <- as.data.frame(colData(spe_not_teen), row.names = rownames(colData(spe_not_teen)))
neuronchat_adult@meta <- as.data.frame(colData(spe_adult), row.names = rownames(colData(spe_adult)))
neuronchat_not_adult@meta <- as.data.frame(colData(spe_not_adult), row.names = rownames(colData(spe_not_adult)))
neuronchat_elderly@meta <- as.data.frame(colData(spe_elderly), row.names = rownames(colData(spe_elderly)))
neuronchat_not_elderly@meta <- as.data.frame(colData(spe_not_elderly), row.names = rownames(colData(spe_not_elderly)))

# Merge NeuronChat objects
object.list <- list(
    Infant = neuronchat_infant,
    Not_Infant = neuronchat_not_infant,
    Teen = neuronchat_teen,
    Not_teen = neuronchat_not_teen,
    Adult = neuronchat_adult,
    Not_adult = neuronchat_not_adult,
    Elderly = neuronchat_elderly,
    Not_elderly = neuronchat_not_elderly
)

neuronchat <-
    mergeNeuronChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Save CellChat objects
save(
    neuronchat_infant,
    neuronchat_not_infant,
    neuronchat_teen,
    neuronchat_not_teen,
    neuronchat_adult,
    neuronchat_not_adult,
    neuronchat_elderly,
    neuronchat_not_elderly,
    neuronchat,
    file = here::here("processed-data", "LR_interactions", "NeuronChat.Rdata")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
