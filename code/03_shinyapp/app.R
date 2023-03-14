library("spatialLIBD")
library("markdown") ## due to a shinyapps.io bug

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data (all paths are relative to this script's location)
spe <- readRDS("QCed_spe.rds")
<<<<<<< HEAD
spe$CellCount <- spe$segmentation_info
=======
spe$CellCount <- spe$NBW
>>>>>>> a15aded873b103f13162e665349087c9e927bd3d
vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "spatial_DG_lifespan, Visium",
    spe_discrete_vars = c(vars[grep("10x_", vars)], "ManualAnnotation"),
    spe_continuous_vars = c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio", "CellCount"),
    default_cluster = "10x_graphclust"
)
