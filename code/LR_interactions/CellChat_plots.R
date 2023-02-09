##################################################
# spatial_DG_lifespan project
# CellChat spot-spot communication pair-wise PLOTS
# Anthony Ramnauth, Feb 08 2023
##################################################

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(CellChat)
    library(patchwork)
})

load(here::here("processed-data", "LR_interactions", "CellChat.Rdata"))

######################################
# Visualization of individual age bins
######################################

# Infant

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interations_number&strength_Infant.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat_infant, title.name = "Number of interactions for Infant")

netVisual_heatmap(cellchat_infant, measure = "weight",
    title.name = "Interaction weights/strength for Infant")

netVisual_heatmap(cellchat, comparison = c(2, 1),
    title.name = "Change in number of interactions for Infant vs. Non-Infant")

netVisual_heatmap(cellchat, comparison = c(2, 1), measure = "weight",
    title.name = "Change in interaction weights/strength for Infant vs. Non-Infant")

groupSize_infant <- as.numeric(table(cellchat_infant@idents))

netVisual_circle(
    cellchat_infant@net$count,
    vertex.weight = groupSize_infant,
    weight.scale = T,
    label.edge = F,
    title.name = "Number of interactions for Infant"
)

netVisual_circle(
    cellchat_infant@net$weight,
    vertex.weight = groupSize_infant,
    weight.scale = T,
    label.edge = F,
    title.name = "Interaction weights/strength for Infant"
)

netVisual_diffInteraction(cellchat, comparison = c(2, 1), weight.scale = T,
    title.name = "Change in number of interactions for Infant vs. Non-Infant")

netVisual_diffInteraction(cellchat, comparison = c(2, 1), weight.scale = T,
    measure = "weight",
    title.name = "Change in interaction weights/strength for Infant vs. Non-Infant")

mat_infant_count <- cellchat_infant@net$count
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_infant_count)) {
    mat2_infant_count <-
        matrix(
            0,
            nrow = nrow(mat_infant_count),
            ncol = ncol(mat_infant_count),
            dimnames = dimnames(mat_infant_count)
        )
    mat2_infant_count[i,] <- mat_infant_count[i,]
    netVisual_circle(
        mat2_infant_count,
        vertex.weight = groupSize_infant,
        weight.scale = T,
        edge.weight.max = max(mat_infant_count),
        title.name = paste0(
            "Number of interactions for cluster #",
            rownames(mat_infant_count)[i]
        )
    )
}

mat_infant <- cellchat_infant@net$weight
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_infant)) {
    mat2_infant <-
        matrix(
            0,
            nrow = nrow(mat_infant),
            ncol = ncol(mat_infant),
            dimnames = dimnames(mat_infant)
        )
    mat2_infant[i,] <- mat_infant[i,]
    netVisual_circle(
        mat2_infant,
        vertex.weight = groupSize_infant,
        weight.scale = T,
        edge.weight.max = max(mat_infant),
        title.name = paste0(
            "Strength of interactions for cluster #",
            rownames(mat_infant)[i]
        )
    )
}

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interations_pathways_Infant.pdf"
    ),
    width = 12,
    height = 10
)

netVisual_bubble(
    cellchat_infant,
    sources.use = c(1:4),
    targets.use = c(1:4),
    remove.isolate = FALSE,
    font.size = 8
)

netVisual_chord_gene(
    cellchat_infant,
    sources.use = c(1:4),
    targets.use = c(1:4),
    slot.name = "netP",
    lab.cex = .4
)

infant_out <-
    netAnalysis_signalingRole_heatmap(
        cellchat_infant,
        pattern = "outgoing",
        font.size = 7,
        title = "Infant"
    )

infant_in <-
    netAnalysis_signalingRole_heatmap(
        cellchat_infant,
        pattern = "incoming",
        font.size = 7,
        title = "Infant"
    )

infant_out + infant_in

rankNet(cellchat, comparison = c(1, 2), mode = "comparison",
    stacked = T, do.stat = TRUE)

rankNet(cellchat, comparison = c(1, 2), mode = "comparison",
    stacked = F, do.stat = TRUE)

dev.off()

# Teen

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interations_number&strength_Teen.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat_teen, title.name = "Number of interactions for Teen")

netVisual_heatmap(cellchat_teen, measure = "weight",
    title.name = "Interaction weights/strength for Teen")

netVisual_heatmap(cellchat, comparison = c(4, 3),
    title.name = "Change in number of interactions for Teen vs. Non-Teen")

netVisual_heatmap(cellchat, comparison = c(4, 3), measure = "weight",
    title.name = "Change in interaction weights/strength for Teen vs. Non-Teen")

groupSize_teen <- as.numeric(table(cellchat_teen@idents))

netVisual_circle(
    cellchat_teen@net$count,
    vertex.weight = groupSize_teen,
    weight.scale = T,
    label.edge = F,
    title.name = "Number of interactions for Teen"
)

netVisual_circle(
    cellchat_teen@net$weight,
    vertex.weight = groupSize_teen,
    weight.scale = T,
    label.edge = F,
    title.name = "Interaction weights/strength for Teen"
)

netVisual_diffInteraction(cellchat, comparison = c(4, 3), weight.scale = T,
    title.name = "Change in number of interactions for Teen vs. Non-Teen")

netVisual_diffInteraction(cellchat, comparison = c(4, 3), weight.scale = T,
    measure = "weight",
    title.name = "Change in interaction weights/strength for Teen vs. Non-Teen")

mat_teen_count <- cellchat_teen@net$count
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_teen_count)) {
    mat2_teen_count <-
        matrix(
            0,
            nrow = nrow(mat_teen_count),
            ncol = ncol(mat_teen_count),
            dimnames = dimnames(mat_teen_count)
        )
    mat2_teen_count[i,] <- mat_teen_count[i,]
    netVisual_circle(
        mat2_teen_count,
        vertex.weight = groupSize_teen,
        weight.scale = T,
        edge.weight.max = max(mat_teen_count),
        title.name = paste0(
            "Number of interactions for cluster #",
            rownames(mat_teen_count)[i]
        )
    )
}

mat_teen <- cellchat_teen@net$weight
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_teen)) {
    mat2_teen <-
        matrix(
            0,
            nrow = nrow(mat_teen),
            ncol = ncol(mat_teen),
            dimnames = dimnames(mat_teen)
        )
    mat2_teen[i,] <- mat_teen[i,]
    netVisual_circle(
        mat2_teen,
        vertex.weight = groupSize_teen,
        weight.scale = T,
        edge.weight.max = max(mat_teen),
        title.name = paste0(
            "Strength of interactions for cluster #",
            rownames(mat_teen)[i]
        )
    )
}

dev.off()

pdf(
    file = here::here("plots", "LR_interactions", "Interations_pathways_Teen.pdf"),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat_teen,
    sources.use = c(1:4),
    targets.use = c(1:4),
    remove.isolate = FALSE,
    font.size = 6
)

netVisual_chord_gene(
    cellchat_teen,
    sources.use = c(1:4),
    targets.use = c(1:4),
    slot.name = "netP",
    lab.cex = .4
)

teen_out <-
    netAnalysis_signalingRole_heatmap(
        cellchat_teen,
        pattern = "outgoing",
        font.size = 5,
        title = "Teen"
    )

teen_in <-
    netAnalysis_signalingRole_heatmap(
        cellchat_teen,
        pattern = "incoming",
        font.size = 5,
        title = "Teen"
    )

teen_out + teen_in

rankNet(cellchat, comparison = c(3, 4), mode = "comparison",
    stacked = T, do.stat = TRUE)

rankNet(cellchat, comparison = c(3, 4), mode = "comparison",
    stacked = F, do.stat = TRUE)

dev.off()

# Adult

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interations_number&strength_Adult.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat_adult, title.name = "Number of interactions for Adult")

netVisual_heatmap(cellchat_adult, measure = "weight",
    title.name = "Interaction weights/strength for Adult")

netVisual_heatmap(cellchat, comparison = c(6, 5),
    title.name = "Change in number of interactions for Adult vs. Non-Adult")

netVisual_heatmap(cellchat, comparison = c(6, 5), measure = "weight",
    title.name = "Change in interaction weights/strength for Adult vs. Non-Adult")

groupSize_adult <- as.numeric(table(cellchat_adult@idents))

netVisual_circle(
    cellchat_adult@net$count,
    vertex.weight = groupSize_adult,
    weight.scale = T,
    label.edge = F,
    title.name = "Number of interactions for Adult"
)

netVisual_circle(
    cellchat_adult@net$weight,
    vertex.weight = groupSize_adult,
    weight.scale = T,
    label.edge = F,
    title.name = "Interaction weights/strength for Adult"
)

netVisual_diffInteraction(cellchat, comparison = c(6, 5), weight.scale = T,
    title.name = "Change in number of interactions for Adult vs. Non-Adult")

netVisual_diffInteraction(cellchat, comparison = c(6, 5), weight.scale = T,
    measure = "weight",
    title.name = "Change in interaction weights/strength for Adult vs. Non-Adult")

mat_adult_count <- cellchat_adult@net$count
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_adult_count)) {
    mat2_adult_count <-
        matrix(
            0,
            nrow = nrow(mat_adult_count),
            ncol = ncol(mat_adult_count),
            dimnames = dimnames(mat_adult_count)
        )
    mat2_adult_count[i,] <- mat_adult_count[i,]
    netVisual_circle(
        mat2_adult_count,
        vertex.weight = groupSize_adult,
        weight.scale = T,
        edge.weight.max = max(mat_adult_count),
        title.name = paste0(
            "Number of interactions for cluster #",
            rownames(mat_adult_count)[i]
        )
    )
}

mat_adult <- cellchat_adult@net$weight
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_adult)) {
    mat2_adult <-
        matrix(
            0,
            nrow = nrow(mat_adult),
            ncol = ncol(mat_adult),
            dimnames = dimnames(mat_adult)
        )
    mat2_adult[i,] <- mat_adult[i,]
    netVisual_circle(
        mat2_adult,
        vertex.weight = groupSize_adult,
        weight.scale = T,
        edge.weight.max = max(mat_adult),
        title.name = paste0(
            "Strength of interactions for cluster #",
            rownames(mat_adult)[i]
        )
    )
}

dev.off()

pdf(
    file = here::here("plots", "LR_interactions", "Interations_pathways_Adult.pdf"),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat_adult,
    sources.use = c(1:4),
    targets.use = c(1:4),
    remove.isolate = FALSE,
    font.size = 6
)

netVisual_chord_gene(
    cellchat_adult,
    sources.use = c(1:4),
    targets.use = c(1:4),
    slot.name = "netP",
    lab.cex = .4
)

adult_out <-
    netAnalysis_signalingRole_heatmap(
        cellchat_adult,
        pattern = "outgoing",
        font.size = 5,
        title = "Adult"
    )

adult_in <-
    netAnalysis_signalingRole_heatmap(
        cellchat_adult,
        pattern = "incoming",
        font.size = 5,
        title = "Adult"
    )

adult_out + adult_in

rankNet(cellchat, comparison = c(5, 6), mode = "comparison",
    stacked = T, do.stat = TRUE)

rankNet(cellchat, comparison = c(5, 6), mode = "comparison",
    stacked = F, do.stat = TRUE)

dev.off()

# Elderly

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interations_number&strength_Elderly.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat_elderly, title.name = "Number of interactions for Elderly")

netVisual_heatmap(cellchat_elderly, measure = "weight",
    title.name = "Interaction weights/strength for Elderly")

netVisual_heatmap(cellchat, comparison = c(8, 7),
    title.name = "Change in number of interactions for Elderly vs. Non-Elderly")

netVisual_heatmap(cellchat, comparison = c(8, 7), measure = "weight",
    title.name = "Change in interaction weights/strength for Elderly vs. Non-Elderly")

groupSize_elderly <- as.numeric(table(cellchat_elderly@idents))

netVisual_circle(
    cellchat_elderly@net$count,
    vertex.weight = groupSize_elderly,
    weight.scale = T,
    label.edge = F,
    title.name = "Number of interactions for Elderly"
)

netVisual_circle(
    cellchat_elderly@net$weight,
    vertex.weight = groupSize_elderly,
    weight.scale = T,
    label.edge = F,
    title.name = "Interaction weights/strength for Elderly"
)

netVisual_diffInteraction(cellchat, comparison = c(8, 7), weight.scale = T,
    title.name = "Change in number of interactions for Elderly vs. Non-Elderly")

netVisual_diffInteraction(cellchat, comparison = c(8, 7), weight.scale = T,
    measure = "weight",
    title.name = "Change in interaction weights/strength for Elderly vs. Non-Elderly")

mat_elderly_count <- cellchat_elderly@net$count
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_elderly_count)) {
    mat2_elderly_count <-
        matrix(
            0,
            nrow = nrow(mat_elderly_count),
            ncol = ncol(mat_elderly_count),
            dimnames = dimnames(mat_elderly_count)
        )
    mat2_elderly_count[i,] <- mat_elderly_count[i,]
    netVisual_circle(
        mat2_elderly_count,
        vertex.weight = groupSize_elderly,
        weight.scale = T,
        edge.weight.max = max(mat_elderly_count),
        title.name = paste0(
            "Number of interactions for cluster #",
            rownames(mat_elderly_count)[i]
        )
    )
}

mat_elderly <- cellchat_elderly@net$weight
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:nrow(mat_elderly)) {
    mat2_elderly <-
        matrix(
            0,
            nrow = nrow(mat_elderly),
            ncol = ncol(mat_elderly),
            dimnames = dimnames(mat_elderly)
        )
    mat2_elderly[i,] <- mat_elderly[i,]
    netVisual_circle(
        mat2_elderly,
        vertex.weight = groupSize_elderly,
        weight.scale = T,
        edge.weight.max = max(mat_elderly),
        title.name = paste0(
            "Strength of interactions for cluster #",
            rownames(mat_elderly)[i]
        )
    )
}

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interations_pathways_Elderly.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat_elderly,
    sources.use = c(1:4),
    targets.use = c(1:4),
    remove.isolate = FALSE,
    font.size = 6
)

netVisual_chord_gene(
    cellchat_elderly,
    sources.use = c(1:4),
    targets.use = c(1:4),
    slot.name = "netP",
    lab.cex = .4
)

elderly_out <-
    netAnalysis_signalingRole_heatmap(
        cellchat_elderly,
        pattern = "outgoing",
        font.size = 5,
        title = "Elderly"
    )

elderly_in <-
    netAnalysis_signalingRole_heatmap(
        cellchat_elderly,
        pattern = "incoming",
        font.size = 5,
        title = "Elderly"
    )

elderly_out + elderly_in

rankNet(cellchat, comparison = c(7, 8), mode = "comparison",
    stacked = T, do.stat = TRUE)

rankNet(cellchat, comparison = c(7, 8), mode = "comparison",
    stacked = F, do.stat = TRUE)

dev.off()

#################################
# Visualization of merged dataset
#################################

## Infant to Teen

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interation_changes_InfanttoTeen.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat, comparison = c(1, 2),
    title.name = "Change in number of interactions from Infant to Teen")

netVisual_heatmap(
    cellchat,
    comparison = c(1, 2),
    measure = "weight",
    title.name = "Change in strength of  interactions from Infant to Teen"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(1, 2),
    weight.scale = T,
    title.name = "Change in number of interactions from Infant to Teen"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(1, 2),
    weight.scale = T,
    measure = "weight",
    title.name = "Change in strength of interactions from Infant to Teen"
)

infant_teen_count <-
    compareInteractions(cellchat, show.legend = F, group = c(1, 2))

infant_teen_strength <-
    compareInteractions(
        cellchat,
        show.legend = F,
        group = c(1, 2),
        measure = "weight"
    )

infant_teen_count + infant_teen_strength

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Pathway_changes_InfanttoTeen.pdf"
    ),
    width = 8,
    height = 10
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(1, 2),
    stacked = T,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Infant to Teen"
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(1, 2),
    stacked = F,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Infant to Teen"
)

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "LR_changes_InfanttoTeen.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat,
    sources.use = c(1:4),
    targets.use = c(1:4),
    comparison = c(1, 2),
    font.size = 6,
    angle.x = 45
)

dev.off()

## Infant to Adult

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interation_changes_InfanttoAdult.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat, comparison = c(1, 3),
    title.name = "Change in number of interactions from Infant to Adult")

netVisual_heatmap(
    cellchat,
    comparison = c(1, 3),
    measure = "weight",
    title.name = "Change in strength of  interactions from Infant to Adult"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(1, 3),
    weight.scale = T,
    title.name = "Change in number of interactions from Infant to Adult"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(1, 3),
    weight.scale = T,
    measure = "weight",
    title.name = "Change in strength of interactions from Infant to Adult"
)

infant_adult_count <-
    compareInteractions(cellchat, show.legend = F, group = c(1, 3))

infant_adult_strength <-
    compareInteractions(
        cellchat,
        show.legend = F,
        group = c(1, 3),
        measure = "weight"
    )

infant_adult_count + infant_adult_strength

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Pathway_changes_InfanttoAdult.pdf"
    ),
    width = 8,
    height = 10
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(1, 3),
    stacked = T,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Infant to Adult"
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(1, 3),
    stacked = F,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Infant to Adult"
)

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "LR_changes_InfanttoAdult.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat,
    sources.use = c(1:4),
    targets.use = c(1:4),
    comparison = c(1, 3),
    font.size = 6,
    angle.x = 45
)

dev.off()

## Infant to Elderly

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interation_changes_InfanttoElderly.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat, comparison = c(1, 4),
    title.name = "Change in number of interactions from Infant to Elderly")

netVisual_heatmap(
    cellchat,
    comparison = c(1, 4),
    measure = "weight",
    title.name = "Change in strength of  interactions from Infant to Elderly"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(1, 4),
    weight.scale = T,
    title.name = "Change in number of interactions from Infant to Elderly"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(1, 4),
    weight.scale = T,
    measure = "weight",
    title.name = "Change in strength of interactions from Infant to Elderly"
)

infant_elderly_count <-
    compareInteractions(cellchat, show.legend = F, group = c(1, 4))

infant_elderly_strength <-
    compareInteractions(
        cellchat,
        show.legend = F,
        group = c(1, 4),
        measure = "weight"
    )

infant_elderly_count + infant_elderly_strength

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Pathway_changes_InfanttoElderly.pdf"
    ),
    width = 8,
    height = 10
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(1, 4),
    stacked = T,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Infant to Elderly"
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(1, 4),
    stacked = F,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Infant to Elderly"
)

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "LR_changes_InfanttoElderly.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat,
    sources.use = c(1:4),
    targets.use = c(1:4),
    comparison = c(1, 4),
    font.size = 6,
    angle.x = 45
)

dev.off()

## Teen to Adult

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interation_changes_TeentoAdult.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat, comparison = c(2, 3),
    title.name = "Change in number of interactions from Teen to Adult")

netVisual_heatmap(
    cellchat,
    comparison = c(2, 3),
    measure = "weight",
    title.name = "Change in strength of  interactions from Teen to Adult"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(2, 3),
    weight.scale = T,
    title.name = "Change in number of interactions from Teen to Adult"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(2, 3),
    weight.scale = T,
    measure = "weight",
    title.name = "Change in strength of interactions from Teen to Adult"
)

teen_adult_count <-
    compareInteractions(cellchat, show.legend = F, group = c(2, 3))

teen_adult_strength <-
    compareInteractions(
        cellchat,
        show.legend = F,
        group = c(2, 3),
        measure = "weight"
    )

teen_adult_count + teen_adult_strength

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Pathway_changes_TeentoAdult.pdf"
    ),
    width = 8,
    height = 10
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(2, 3),
    stacked = T,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Teen to Adult"
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(2, 3),
    stacked = F,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Teen to Adult"
)

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "LR_changes_TeentoAdult.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat,
    sources.use = c(1:4),
    targets.use = c(1:4),
    comparison = c(2, 3),
    font.size = 6,
    angle.x = 45
)

dev.off()

## Teen to Elderly

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interation_changes_TeentoElderly.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat, comparison = c(2, 4),
    title.name = "Change in number of interactions from Teen to Elderly")

netVisual_heatmap(
    cellchat,
    comparison = c(2, 4),
    measure = "weight",
    title.name = "Change in strength of  interactions from Teen to Elderly"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(2, 4),
    weight.scale = T,
    title.name = "Change in number of interactions from Teen to Elderly"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(2, 4),
    weight.scale = T,
    measure = "weight",
    title.name = "Change in strength of interactions from Teen to Elderly"
)

teen_elderly_count <-
    compareInteractions(cellchat, show.legend = F, group = c(2, 4))

teen_elderly_strength <-
    compareInteractions(
        cellchat,
        show.legend = F,
        group = c(2, 4),
        measure = "weight"
    )

teen_elderly_count + teen_elderly_strength

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Pathway_changes_TeentoElderly.pdf"
    ),
    width = 8,
    height = 10
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(2, 4),
    stacked = T,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Teen to Elderly"
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(2, 4),
    stacked = F,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Teen to Elderly"
)

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "LR_changes_TeentoElderly.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat,
    sources.use = c(1:4),
    targets.use = c(1:4),
    comparison = c(2, 4),
    font.size = 6,
    angle.x = 45
)

dev.off()

## Adult to Elderly

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Interation_changes_AdulttoElderly.pdf"
    ),
    width = 8,
    height = 6
)

netVisual_heatmap(cellchat, comparison = c(3, 4),
    title.name = "Change in number of interactions from Adult to Elderly")

netVisual_heatmap(
    cellchat,
    comparison = c(3, 4),
    measure = "weight",
    title.name = "Change in strength of  interactions from Adult to Elderly"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(3, 4),
    weight.scale = T,
    title.name = "Change in number of interactions from Adult to Elderly"
)

netVisual_diffInteraction(
    cellchat,
    comparison = c(3, 4),
    weight.scale = T,
    measure = "weight",
    title.name = "Change in strength of interactions from Adult to Elderly"
)

adult_elderly_count <-
    compareInteractions(cellchat, show.legend = F, group = c(3, 4))

adult_elderly_strength <-
    compareInteractions(
        cellchat,
        show.legend = F,
        group = c(3, 4),
        measure = "weight"
    )

adult_elderly_count + adult_elderly_strength

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Pathway_changes_AdulttoElderly.pdf"
    ),
    width = 8,
    height = 10
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(3, 4),
    stacked = T,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Adult to Elderly"
)

rankNet(
    cellchat,
    mode = "comparison",
    comparison = c(3, 4),
    stacked = F,
    do.stat = TRUE,
    title = "Changed or conserved pathways from Adult to Elderly"
)

dev.off()

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "LR_changes_AdulttoElderly.pdf"
    ),
    width = 12,
    height = 12
)

netVisual_bubble(
    cellchat,
    sources.use = c(1:4),
    targets.use = c(1:4),
    comparison = c(3, 4),
    font.size = 6,
    angle.x = 45
)

dev.off()



## Reproducibility informationSE)
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
