###################################################
# spatial_DG_lifespan project
# nnSVG per Capture Area
# Anthony Ramnauth, June 29 2022
###################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(scran)
    library(nnSVG)
    library(ggplot2)
    library(here)
    library(sessioninfo)
})

spe <- readRDS(here::here("processed-data", "QC_processed_spe", "QCed_spe.rds"))

## subset spe data based on Capture Area
Br1412_spe <- spe[, spe$sample_id %in% c("Br1412")]

# filter low-expressed and mitochondrial genes
# using default filtering parameters
Br1412_spe <- filter_genes(Br1412_spe)

# set seed for reproducibility
set.seed(12345)

# using a single thread in this example
Br1412_spe <- nnSVG(Br1412_spe)

# number of significant SVGs
table(rowData(Br1412_spe)$padj <= 0.05)

# show results for top n SVGs
rowData(Br1412_spe)[order(rowData(Br1412_spe)$rank)[1:10], ]

# plot spatial expression of top-ranked SVG
ix <- which(rowData(Br1412_spe)$rank == 1)

ix_name <- rowData(Br1412_spe)$gene_name[ix]

ix_name

df <- as.data.frame(
  cbind(spatialCoords(Br1412_spe),
        expr = counts(Br1412_spe)[ix, ]))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = expr)) +
  geom_point(size = 0.8) +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_gradient(low = "gray90", high = "blue", name = "counts") +
  ggtitle(ix_name) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
