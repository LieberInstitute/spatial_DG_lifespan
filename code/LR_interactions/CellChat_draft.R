######################################
# spatial_DG_lifespan project
# CellChat for spot-spot communication
# Anthony Ramnauth, Feb 02 2023
######################################

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
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "clustering_results"),
    prefix = ""
)

# Subset for single capture area
spe_Br1412 <- spe[,spe$sample_id == "Br1412"]

# Set rownames as gene names since ensemble IDs seem not to work
rownames(spe_Br1412) <- rowData(spe_Br1412)$gene_name

# Add imaging scale factors
Br1412_scale.factors <- jsonlite::fromJSON(txt = here::here("processed-data", "01_spaceranger", "Br1412",
    "outs", "spatial", "scalefactors_json.json"))

Br1412_scale.factors <- list(spot.diameter = 55,
    spot = Br1412_scale.factors$spot_diameter_fullres, # these two information are required
    fiducial = Br1412_scale.factors$fiducial_diameter_fullres,
    hires = Br1412_scale.factors$tissue_hires_scalef,
    lowres = Br1412_scale.factors$tissue_lowres_scalef # these three information are not required
)

# Create cellchat object
cellchat <- createCellChat(object = spe_Br1412,
    group.by = "bayesSpace_harmony_8",
    datatype = "spatial",
    coordinates = spatialCoords(spe_Br1412),
    scale.factors = Br1412_scale.factors
)

# Set human database for use
CellChatDB <- CellChatDB.human

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Add to cellchat object
cellchat@DB <- CellChatDB.use

# Preprocess

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in
#order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

# Infer intercellular communication network of each ligand-receptor pair
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                               distance.use = TRUE, interaction.length = 200, scale.distance = 0.01, raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F,
    title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F,
    title.name = "Interaction weights/strength")

# Visualize interactions & strengths
mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Show one pathway
pathways.show <- c("OPIOID")
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)

#visualize the cell-cell communication mediated by a single ligand-receptor pair
pairLR.OPIOID <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)

LR.show <- pairLR.OPIOID[1,]

# Hierarchy plot
vertex.receiver = seq(1,4)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use')
# to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,4,8), remove.isolate = FALSE) # Too much info

# show all the significant unidirectional interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,4,8), signaling = c("OPIOID"), remove.isolate = FALSE)

# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "OPIOID")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the cell-cell communication networks of interest
netAnalysis_signalingRole_scatter(cellchat, signaling = c("OPIOID"))

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
gg2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
gg1+gg2

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
