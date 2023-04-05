######################################
# spatial_DG_lifespan project
# PRECAST k5 to k15
# Anthony Ramnauth, March 30 2023
######################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
  library("dplyr")
  library("purrr")
  library("Seurat")
  library("here")
  library("sessioninfo")
  library("SpatialExperiment")
  library("PRECAST")
  library("tictoc")
})

load(file = here::here("processed-data", "Seurat", "seuList.Rdata"))

preobj = CreatePRECASTObject(seuList = seuList, selectGenesMethod="HVGs")
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 8, maxIter = 30, verbose = TRUE)

K <- as.numeric(Sys.getenv("SGE_TASK_ID"))

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

save(PRECASTObj, file = here::here("processed-data", "Seurat", paste0("PRECASTObj_",K,".Rdata")))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
