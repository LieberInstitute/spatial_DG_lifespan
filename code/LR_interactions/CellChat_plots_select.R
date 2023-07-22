###########################################################
# spatial_DG_lifespan project
# Select CellChat & NeuronChat plots fpr prolif & inflamm
# Anthony Ramnauth, July 06 2023
###########################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(CellChat)
    library(NeuronChat)
    library(patchwork)
})

load(here::here("processed-data", "LR_interactions", "CellChat.Rdata"))

load(here::here("processed-data", "LR_interactions", "NeuronChat.Rdata"))

neuronchat_DG <-
    readRDS(here::here("processed-data", "LR_interactions", "NeuronChat_DG.rds"))

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Info_flow_InfantEnriched_comparison.pdf"
    ))

infant_signals <- c("GAS", "PERIOSTIN", "PDGF", "AGRN", "PTN", "NOTCH")

rankNet(cellchat, comparison = c(1, 2), mode = "comparison",
    signaling = infant_signals, stacked = T, do.stat = TRUE, font.size = 20) +
    theme(text = element_text(size = 20))

rankNet(cellchat, comparison = c(3, 4), mode = "comparison",
    signaling = infant_signals, stacked = T, do.stat = TRUE, font.size = 20) +
    theme(text = element_text(size = 20))

rankNet(cellchat, comparison = c(5, 6), mode = "comparison",
    signaling = infant_signals, stacked = T, do.stat = TRUE, font.size = 20) +
    theme(text = element_text(size = 20))

rankNet(cellchat, comparison = c(7, 8), mode = "comparison",
    signaling = infant_signals, stacked = T, do.stat = TRUE, font.size = 20) +
    theme(text = element_text(size = 20))

netVisual_bubble(
    cellchat_infant,
    signaling = infant_signals,
    sources.use = c(1:4),
    targets.use = c(1:4),
    remove.isolate = FALSE,
    angle.x = 45,
    font.size = 14
)

dev.off()

######################################################################################

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Info_flow_Neuron_SST.pdf"
    ))

sst <- c("SST_SSTR1", "SST_SSTR2", "SST_SSTR3", "SST_SSTR4")

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(1,2),
    pairLR = sst, stacked = TRUE, do.stat = FALSE, font.size = 20) +
    theme(text = element_text(size = 20))

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(3,4),
    pairLR = sst, stacked = TRUE, do.stat = FALSE, font.size = 20) +
    theme(text = element_text(size = 20))

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(5,6),
    pairLR = sst, stacked = TRUE, do.stat = FALSE, font.size = 20) +
    theme(text = element_text(size = 20))

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(7,8),
    pairLR = sst, stacked = TRUE, do.stat = FALSE, font.size = 20) +
    theme(text = element_text(size = 20))

dev.off()

