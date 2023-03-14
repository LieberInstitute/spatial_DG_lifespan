setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

library("spatialLIBD")
library("lobstr")
library("here")

# Load SPE
spe <- readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

## Check how big it is in memory
lobstr::obj_size(spe)
# 5.22 GB
## That's too big for shinyapps.io. Aim to have an object near 2GB.

## Subset the spe object outside of shinyapps.io. Otherwise, the peak memory is
## still affected by loading the object.
## Also, running lobstr::obj_size() takes a while to run, which we don't need
## to run every time someone accesses the shiny app.
imgData(spe) <-
    imgData(spe)[!imgData(spe)$image_id %in% c("hires", "detected", "aligned"), ]
assays(spe)$counts <- NULL
lobstr::obj_size(spe)
# 1.86 GB
## Ok, this is reasonable.

## Save the reduced version of the spe object in the shiny app directory
## instead of using soft links.
saveRDS(spe, file = here::here("code", "03_shinyapp", "QCed_spe.rds"))
