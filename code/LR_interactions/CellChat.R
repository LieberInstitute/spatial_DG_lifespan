#########################################
# spatial_DG_lifespan project
# CellChat spot-spot communication
# Anthony Ramnauth, Feb 02 2023
#########################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(SpatialExperiment)
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
spe_infant <- spe_infant[, spe_infant$bayesSpace_harmony_10 == "2" |
        spe_infant$bayesSpace_harmony_10 == "4" |
        spe_infant$bayesSpace_harmony_10 == "6" |
        spe_infant$bayesSpace_harmony_10 == "87"]

spe_not_infant <- spe_not_infant[, spe_not_infant$bayesSpace_harmony_10 == "2" |
        spe_not_infant$bayesSpace_harmony_10 == "4" |
        spe_not_infant$bayesSpace_harmony_10 == "6" |
        spe_not_infant$bayesSpace_harmony_10 == "7"]

spe_teen <- spe_teen[, spe_teen$bayesSpace_harmony_10 == "2" |
        spe_teen$bayesSpace_harmony_10 == "4" |
        spe_teen$bayesSpace_harmony_10 == "6" |
        spe_teen$bayesSpace_harmony_10 == "7"]

spe_not_teen <- spe_not_teen[, spe_not_teen$bayesSpace_harmony_10 == "2" |
        spe_not_teen$bayesSpace_harmony_10 == "4" |
        spe_not_teen$bayesSpace_harmony_10 == "6" |
        spe_not_teen$bayesSpace_harmony_10 == "7"]

spe_adult <- spe_adult[, spe_adult$bayesSpace_harmony_10 == "2" |
        spe_adult$bayesSpace_harmony_10 == "4" |
        spe_adult$bayesSpace_harmony_10 == "6" |
        spe_adult$bayesSpace_harmony_10 == "7"]

spe_not_adult <- spe_not_adult[, spe_not_adult$bayesSpace_harmony_10 == "2" |
        spe_not_adult$bayesSpace_harmony_10 == "4" |
        spe_not_adult$bayesSpace_harmony_10 == "6" |
        spe_not_adult$bayesSpace_harmony_10 == "7"]

spe_elderly <- spe_elderly[, spe_elderly$bayesSpace_harmony_10 == "2" |
            spe_elderly$bayesSpace_harmony_10 == "4" |
            spe_elderly$bayesSpace_harmony_10 == "6" |
            spe_elderly$bayesSpace_harmony_10 == "7"]

spe_not_elderly <- spe_not_elderly[, spe_not_elderly$bayesSpace_harmony_10 == "2" |
            spe_not_elderly$bayesSpace_harmony_10 == "4" |
            spe_not_elderly$bayesSpace_harmony_10 == "6" |
            spe_not_elderly$bayesSpace_harmony_10 == "7"]

# Create cellchat object for age binned spe objects
cellchat_infant <- createCellChat(object = spe_infant,
    group.by = "BayesSpace")

cellchat_not_infant <- createCellChat(object = spe_not_infant,
    group.by = "BayesSpace")

cellchat_teen <- createCellChat(object = spe_teen,
    group.by = "BayesSpace")

cellchat_not_teen <- createCellChat(object = spe_not_teen,
    group.by = "BayesSpace")

cellchat_adult <- createCellChat(object = spe_adult,
    group.by = "BayesSpace")

cellchat_not_adult <- createCellChat(object = spe_not_adult,
    group.by = "BayesSpace")

cellchat_elderly <- createCellChat(object = spe_elderly,
    group.by = "BayesSpace")

cellchat_not_elderly <- createCellChat(object = spe_not_elderly,
    group.by = "BayesSpace")

# Set human database for use
CellChatDB <- CellChatDB.human

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Add to cellchat object
cellchat_infant@DB <- CellChatDB.use
cellchat_not_infant@DB <- CellChatDB.use
cellchat_teen@DB <- CellChatDB.use
cellchat_not_teen@DB <- CellChatDB.use
cellchat_adult@DB <- CellChatDB.use
cellchat_not_adult@DB <- CellChatDB.use
cellchat_elderly@DB <- CellChatDB.use
cellchat_not_elderly@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_infant <- subsetData(cellchat_infant)
cellchat_not_infant <- subsetData(cellchat_not_infant)
cellchat_teen <- subsetData(cellchat_teen)
cellchat_not_teen <- subsetData(cellchat_not_teen)
cellchat_adult <- subsetData(cellchat_adult)
cellchat_not_adult <- subsetData(cellchat_not_adult)
cellchat_elderly <- subsetData(cellchat_elderly)
cellchat_not_elderly <- subsetData(cellchat_not_elderly)

#future::plan("multiprocess", workers = 4) # do parallel (optional?)

cellchat_infant <- identifyOverExpressedGenes(cellchat_infant)
cellchat_not_infant <- identifyOverExpressedGenes(cellchat_not_infant)
cellchat_teen <- identifyOverExpressedGenes(cellchat_teen)
cellchat_not_teen <- identifyOverExpressedGenes(cellchat_not_teen)
cellchat_adult <- identifyOverExpressedGenes(cellchat_adult)
cellchat_not_adult <- identifyOverExpressedGenes(cellchat_not_adult)
cellchat_elderly <- identifyOverExpressedGenes(cellchat_elderly)
cellchat_not_elderly <- identifyOverExpressedGenes(cellchat_not_elderly)

cellchat_infant <- identifyOverExpressedInteractions(cellchat_infant)
cellchat_not_infant <- identifyOverExpressedInteractions(cellchat_not_infant)
cellchat_teen <- identifyOverExpressedInteractions(cellchat_teen)
cellchat_not_teen <- identifyOverExpressedInteractions(cellchat_not_teen)
cellchat_adult <- identifyOverExpressedInteractions(cellchat_adult)
cellchat_not_adult <- identifyOverExpressedInteractions(cellchat_not_adult)
cellchat_elderly <- identifyOverExpressedInteractions(cellchat_elderly)
cellchat_not_elderly <- identifyOverExpressedInteractions(cellchat_not_elderly)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE`
#in the function `computeCommunProb()` in
#order to use the projected data)
cellchat_infant <- projectData(cellchat_infant, PPI.human)
cellchat_not_infant <- projectData(cellchat_not_infant, PPI.human)
cellchat_teen <- projectData(cellchat_teen, PPI.human)
cellchat_not_teen <- projectData(cellchat_not_teen, PPI.human)
cellchat_adult <- projectData(cellchat_adult, PPI.human)
cellchat_not_adult <- projectData(cellchat_not_adult, PPI.human)
cellchat_elderly <- projectData(cellchat_elderly, PPI.human)
cellchat_not_elderly <- projectData(cellchat_not_elderly, PPI.human)

# Infer intercellular communication network of each ligand-receptor pair
cellchat_infant <-
    computeCommunProb(cellchat_infant,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_not_infant <-
    computeCommunProb(cellchat_not_infant,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_teen <-
    computeCommunProb(cellchat_teen,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_not_teen <-
    computeCommunProb(cellchat_not_teen,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_adult <-
    computeCommunProb(cellchat_adult,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_not_adult <-
    computeCommunProb(cellchat_not_adult,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_elderly <-
    computeCommunProb(cellchat_elderly,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

cellchat_not_elderly <-
    computeCommunProb(cellchat_not_elderly,
        type = "truncatedMean",
        trim = 0.1,
        raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_infant <- filterCommunication(cellchat_infant, min.cells = 10)
cellchat_not_infant <- filterCommunication(cellchat_not_infant, min.cells = 10)
cellchat_teen <- filterCommunication(cellchat_teen, min.cells = 10)
cellchat_not_teen <- filterCommunication(cellchat_not_teen, min.cells = 10)
cellchat_adult <- filterCommunication(cellchat_adult, min.cells = 10)
cellchat_not_adult <- filterCommunication(cellchat_not_adult, min.cells = 10)
cellchat_elderly <- filterCommunication(cellchat_elderly, min.cells = 10)
cellchat_not_elderly <- filterCommunication(cellchat_not_elderly, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat_infant <- computeCommunProbPathway(cellchat_infant)
cellchat_not_infant <- computeCommunProbPathway(cellchat_not_infant)
cellchat_teen <- computeCommunProbPathway(cellchat_teen)
cellchat_not_teen <- computeCommunProbPathway(cellchat_not_teen)
cellchat_adult <- computeCommunProbPathway(cellchat_adult)
cellchat_not_adult <- computeCommunProbPathway(cellchat_not_adult)
cellchat_elderly <- computeCommunProbPathway(cellchat_elderly)
cellchat_not_elderly <- computeCommunProbPathway(cellchat_not_elderly)

# Calculate the aggregated cell-cell communication network
cellchat_infant <- aggregateNet(cellchat_infant)
cellchat_not_infant <- aggregateNet(cellchat_not_infant)
cellchat_teen <- aggregateNet(cellchat_teen)
cellchat_not_teen <- aggregateNet(cellchat_not_teen)
cellchat_adult <- aggregateNet(cellchat_adult)
cellchat_not_adult <- aggregateNet(cellchat_not_adult)
cellchat_elderly <- aggregateNet(cellchat_elderly)
cellchat_not_elderly <- aggregateNet(cellchat_not_elderly)

# Compute the network centrality scores
cellchat_infant <- netAnalysis_computeCentrality(cellchat_infant, slot.name = "netP")
cellchat_not_infant <- netAnalysis_computeCentrality(cellchat_not_infant, slot.name = "netP")
cellchat_teen <- netAnalysis_computeCentrality(cellchat_teen, slot.name = "netP")
cellchat_not_teen <- netAnalysis_computeCentrality(cellchat_not_teen, slot.name = "netP")
cellchat_adult <- netAnalysis_computeCentrality(cellchat_adult, slot.name = "netP")
cellchat_not_adult <- netAnalysis_computeCentrality(cellchat_not_adult, slot.name = "netP")
cellchat_elderly <- netAnalysis_computeCentrality(cellchat_elderly, slot.name = "netP")
cellchat_elderly <- netAnalysis_computeCentrality(cellchat_elderly, slot.name = "netP")

################
# Merge datasets
################

# Merge CellChat objects
object.list <- list(
    Infant = cellchat_infant,
    Not_Infant = cellchat_not_infant,
    Teen = cellchat_teen,
    Not_teen = cellchat_not_teen,
    Adult = cellchat_adult,
    Not_adult = cellchat_not_adult,
    Elderly = cellchat_elderly,
    Not_elderly = cellchat_not_elderly
)

cellchat <-
    mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Save CellChat objects
save(
    cellchat_infant,
    cellchat_not_infant,
    cellchat_teen,
    cellchat_not_teen,
    cellchat_adult,
    cellchat_not_adult,
    cellchat_elderly,
    cellchat_not_elderly,
    cellchat,
    file = here::here("processed-data", "LR_interactions", "CellChat.Rdata")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
