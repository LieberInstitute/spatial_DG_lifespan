#############################
# Spatial_DG_lifespan project
# nnSVG results
# Anthony Ramnauth, July 2022
#############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

# load spe object
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# load results from previous script
res_list <- readRDS(here::here("processed-data", "nnSVG", "whole_tissue", "DG_nnSVG_results.rds"))

names(res_list)

# Number of significant SVGs per capture area

table(res_list$Br1412$padj <= 0.05)
# FALSE  TRUE
#  275  4404

table(res_list$Br2706$padj <= 0.05)
# FALSE  TRUE
#   111  1358

table(res_list$Br3942$padj <= 0.05)
# FALSE  TRUE
#   189  1589

table(res_list$Br5242$padj <= 0.05)
# FALSE  TRUE
#   428  2239

table(res_list$Br6023$padj <= 0.05)
# FALSE  TRUE
#   297  1997

table(res_list$Br8195$padj <= 0.05)
# FALSE  TRUE
#   234  2068

table(res_list$Br8667$padj <= 0.05)
# FALSE  TRUE
#   154  1956

table(res_list$Br8686$padj <= 0.05)
# FALSE  TRUE
#   165   830


# Create vector of samples for nnSVG on whole tissue
sample_ids <- c(
    "Br1412",
    "Br2706",
    "Br3942",
    "Br5242",
    "Br6023",
    "Br8195",
    "Br8667",
    "Br8686"
)

# ---------------
# combine results
# ---------------

# sum gene ranks across sample-parts to generate overall ranking

# number of genes that passed filtering for each sample-part
sapply(res_list, nrow)

# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_ids

for (s in seq_along(sample_ids)) {
  stopifnot(colnames(res_ranks)[s] == sample_ids[s])
  stopifnot(colnames(res_ranks)[s] == names(res_list)[s])

  rownames_s <- rownames(res_list[[s]])
  res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}

# keep only genes that were not filtered out in all samples
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))

# summary table
df_summary <- data.frame(
  gene_id = names(avg_ranks),
  gene_name = rowData(spe)[names(avg_ranks), "gene_name"],
  gene_type = rowData(spe)[names(avg_ranks), "gene_type"],
  avg_rank = unname(avg_ranks),
  row.names = names(avg_ranks)
)

head(df_summary, 20)

# directory to save whole tissue results
dir_outputs <- here("processed-data", "nnSVG", "whole_tissue")
fn_out <- file.path(dir_outputs, "DG_nnSVG_avgrank")

# Export summary as .csv file
write.csv(df_summary,fn_out, row.names = FALSE)

# Plot top 20 SVGs
SVGs <- c(
    "MBP",
    "PLP1",
    "GFAP",
    "MTRNR2L12",
    "TTR",
    "SNAP25",
    "SLC17A7",
    "UCHL1",
    "NPTXR",
    "THY1",
    "HPCA",
    "ENC1",
    "NPTX1",
    "CHN1",
    "NRGN",
    "YWHAH",
    "OLFM1",
    "NNAT",
    "NCDN",
    "CRYAB"
)

# Locate the marker genes
SVG_search <- rowData(spe)$gene_search[match(SVGs, rowData(spe)$gene_name)]

for (i in SVG_search) {
  vis_grid_gene(
    spe = spe,
    geneid = i,
    pdf = here::here("plots", "nnSVG", "whole_tissue", paste0(gsub("; ", "_", i), ".pdf")),
    assayname = "logcounts",
    minCount = 0,
    cont_colors = c("aquamarine4",
    "springgreen", "goldenrod", "red"),
    alpha = 0.5,
    sample_order = unique(spe$sample_id),
    point_size = 2
  )
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
