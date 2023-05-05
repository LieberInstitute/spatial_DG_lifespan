##################################################
# spatial_DG_lifespan project
# NeuronChat spot-spot communication PLOTS
# Anthony Ramnauth, May 05 2023
##################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(NeuronChat)
    library(CellChat)
    library(patchwork)
})

load(here::here("processed-data", "LR_interactions", "NeuronChat.Rdata"))

neuronchat_DG <-
    readRDS(here::here("processed-data", "LR_interactions", "NeuronChat_DG.rds"))

#######
#all DG
#######


pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Neuron_LR_All.pdf"
    ),
    width = 14,
    height = 10
)


rankNet_Neuron(neuronchat_DG, mode = "single", measure = c("weight"),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat_DG, mode = "single", measure = c("weight"),
    stacked = FALSE ,do.stat = TRUE)

rankNet_Neuron(neuronchat_DG, mode = "single", measure = c("count"),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat_DG, mode = "single", measure = c("count"),
    stacked = FALSE ,do.stat = TRUE)

lig_tar_heatmap(neuronchat_DG, interaction_name = "GJA1_GJA1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIA2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIA1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIN2A", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIA3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIN2B", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIK5", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN1_NLGN2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN1_NLGN3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN2_NLGN2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN1_NLGN4X", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRM5", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN2_NLGN3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN3_NLGN2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN3_NLGN3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRM3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIA4", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN2_NLGN4X", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN3_NLGN4X", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRM7", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "CO_GUCY1B1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN1_NLGN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "SST_SSTR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRM2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIK2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN2_NLGN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRM1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "Glu_GRIK4", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NRXN3_NLGN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NPY_NPY1R", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "SST_SSTR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "CO_GUCY1A2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "NPY_NPY2R", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "PDYN_OPRD1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "CORT_SSTR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "GJB6_GJB6", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_DG, interaction_name = "CORT_SSTR1", width.vector=c(0.38,0.35,0.27))


dev.off()


######################################
# Visualization of individual age bins
######################################

# Infant


pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Neuron_LR_Infant.pdf"
    ),
    width = 14,
    height = 10
)


rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(1,2),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(1,2),
    stacked = FALSE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(1,2),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(1,2),
    stacked = FALSE ,do.stat = TRUE)

lig_tar_heatmap(neuronchat_infant, interaction_name = "Glu_GRIN1", width.vector=c(0.38,0.35,0.27))

dev.off()

# Teen

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Neuron_LR_Teen.pdf"
    ),
    width = 14,
    height = 10
)


rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(3,4),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(3,4),
    stacked = FALSE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(3,4),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(3,4),
    stacked = FALSE ,do.stat = TRUE)

lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABBR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRA5", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRG2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRB3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABBR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRA1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRB2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRA2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "PDYN_OPRD1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRG1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "Glu_GRIN2C", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "CORT_SSTR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "GABA_GABRA4", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "CO_GUCY1A1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "CORT_SSTR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "SST_SSTR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "NPY_NPY2R", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "NRXN3_NLGN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "NRXN2_NLGN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "NRXN1_NLGN1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "SST_SSTR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "Glu_GRM1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "NRXN3_NLGN2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "CO_GUCY1B1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "Glu_GRM5", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "NPY_NPY1R", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_teen, interaction_name = "Glu_GRIA4", width.vector=c(0.38,0.35,0.27))

dev.off()

# Adult

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Neuron_LR_Adult.pdf"
    ),
    width = 14,
    height = 10
)


rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(5,6),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(5,6),
    stacked = FALSE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(5,6),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(5,6),
    stacked = FALSE ,do.stat = TRUE)

lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABBR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRA5", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRB3", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABBR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRG2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRA1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRB2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRA2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "Glu_GRIN3A", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRG1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "CORT_SSTR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "PDYN_OPRK1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GABA_GABRA4", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "PDYN_OPRD1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "CORT_SSTR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "CO_GUCY1A1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "GJB6_GJB6", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "NPY_NPY1R", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_adult, interaction_name = "Glu_GRIA4", width.vector=c(0.38,0.35,0.27))

dev.off()

# Elderly

pdf(
    file = here::here(
        "plots",
        "LR_interactions",
        "Neuron_LR_Elderly.pdf"
    ),
    width = 14,
    height = 10
)


rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(7,8),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("weight"), comparison = c(7,8),
    stacked = FALSE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(7,8),
    stacked = TRUE ,do.stat = TRUE)

rankNet_Neuron(neuronchat, mode = "comparison", measure = c("count"), comparison = c(7,8),
    stacked = FALSE ,do.stat = TRUE)

lig_tar_heatmap(neuronchat_elderly, interaction_name = "NRXN1_NLGN4Y", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "NRXN2_NLGN4Y", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "NRXN3_NLGN4Y", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "CORT_SSTR2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "PDYN_OPRK1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "PDYN_OPRD1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "CORT_SSTR1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "GJA1_GJA1", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "NRXN3_NLGN2", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "Glu_GRIA4", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "Glu_GRM5", width.vector=c(0.38,0.35,0.27))
lig_tar_heatmap(neuronchat_elderly, interaction_name = "NRXN3_NLGN3", width.vector=c(0.38,0.35,0.27))

dev.off()

## Reproducibility informationSE)
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
