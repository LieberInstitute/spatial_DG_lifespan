################################################
# Spatial_DG_lifespan project
# Finding overlapping markers from various tests
# Anthony Ramnauth, July 19 2022
################################################

library(dplyr)
library(sessioninfo)
library(here)

SVG <- read.csv(here::here("processed-data", "nnSVG", "BayesSpace", "DG_BayesSpace_nnSVG_avgrank.csv"))
binom1 <- read.csv(here::here("processed-data", "BayesSpace", "Clust_1_binomial_test_results.csv"))
enriched1 <- read.csv(here::here("processed-data", "BayesSpace", "Clust_1_enriched_results.csv"))
binom2 <- read.csv(here::here("processed-data", "BayesSpace", "Clust_2_binomial_test_results.csv"))
enriched2 <- read.csv(here::here("processed-data", "BayesSpace", "Clust_2_enriched_results.csv"))
binomGCL <- read.csv(here::here("processed-data", "BayesSpace", "GCL_binomial_test_results.csv"))
enrichedGCL <- read.csv(here::here("processed-data", "BayesSpace", "GCL_enriched_results.csv"))
binomSGZ <- read.csv(here::here("processed-data", "BayesSpace", "SGZ_binomial_test_results.csv"))
enrichedSGZ <- read.csv(here::here("processed-data", "BayesSpace", "SGZ_enriched_results.csv"))
binomCA4 <- read.csv(here::here("processed-data", "BayesSpace", "CA4_binomial_test_results.csv"))
enrichedCA4 <- read.csv(here::here("processed-data", "BayesSpace", "CA4_enriched_results.csv"))
binomCA3 <- read.csv(here::here("processed-data", "BayesSpace", "CA3_binomial_test_results.csv"))
enrichedCA3 <- read.csv(here::here("processed-data", "BayesSpace", "CA3_enriched_results.csv"))
binomML <- read.csv(here::here("processed-data", "BayesSpace", "ML_binomial_test_results.csv"))
enrichedML <- read.csv(here::here("processed-data", "BayesSpace", "ML_enriched_results.csv"))
binom8 <- read.csv(here::here("processed-data", "BayesSpace", "Clust_8_binomial_test_results.csv"))
enriched8 <- read.csv(here::here("processed-data", "BayesSpace", "Clust_8_enriched_results.csv"))

# Create vectors of the gene names
en1 <- c(enriched1$gene_name)
bi1 <- c(binom1$gene_name)
svg_list <- c(SVG$gene_name)

# Find gene names from one vector in the other
bien1 <- en1[en1 %in% bi1]
svg1 <- bien1[bien1 %in% svg_list]
svg1 <- svg1[!is.na(svg1)]

# Create vectors of the gene names
en2 <- c(enriched2$gene_name)
bi2 <- c(binom2$gene_name)

# Find gene names from one vector in the other
bien2 <- en2[en2 %in% bi2]
svg2 <- bien2[bien2 %in% svg_list]
svg2 <- svg2[!is.na(svg2)]

# Create vectors of the gene names
enGCL <- c(enrichedGCL$gene_name)
biGCL <- c(binomGCL$gene_name)

# Find gene names from one vector in the other
bienGCL <- enGCL[enGCL %in% biGCL]
svgGCL <- bienGCL[bienGCL %in% svg_list]
svgGCL <- svgGCL[!is.na(svgGCL)]

# Create vectors of the gene names
enSGZ <- c(enrichedSGZ$gene_name)
biSGZ <- c(binomSGZ$gene_name)

# Find gene names from one vector in the other
bienSGZ <- enSGZ[enSGZ %in% biSGZ]
svgSGZ <- bienSGZ[bienSGZ %in% svg_list]
svgSGZ <- svgSGZ[!is.na(svgSGZ)]

# Create vectors of the gene names
enCA4 <- c(enrichedCA4$gene_name)
biCA4 <- c(binomCA4$gene_name)

# Find gene names from one vector in the other
bienCA4 <- enCA4[enCA4 %in% biCA4]
svgCA4 <- bienCA4[bienCA4 %in% svg_list]
svgCA4 <- svgCA4[!is.na(svgCA4)]

# Create vectors of the gene names
enCA3 <- c(enrichedCA3$gene_name)
biCA3 <- c(binomCA3$gene_name)

# Find gene names from one vector in the other
bienCA3 <- enCA3[enCA3 %in% biCA3]
svgCA3 <- bienCA3[bienCA3 %in% svg_list]
svgCA3 <- svgCA3[!is.na(svgCA3)]

# Create vectors of the gene names
enML <- c(enrichedML$gene_name)
biML <- c(binomML$gene_name)

# Find gene names from one vector in the other
bienML <- enML[enML %in% biML]
svgML <- bienML[bienML %in% svg_list]
svgML <- svgML[!is.na(svgML)]

# Create vectors of the gene names
en8 <- c(enriched8$gene_name)
bi8 <- c(binom8$gene_name)

# Find gene names from one vector in the other
bien8 <- en8[en8 %in% bi8]
svg8 <- bien8[bien8 %in% svg_list]
svg8 <- svg8[!is.na(svg8)]

# Make a dataframe to export to .csv file

max_ln <- max(c(
    length(svg1),
    length(svg2),
    length(svgGCL),
    length(svgSGZ),
    length(svgCA4),
    length(svgCA3),
    length(svgML),
    length(svg8)
    ))

gene_markers <- data.frame(
    clust_1 = c(svg1, rep(NA, max_ln - length(svg1))),
    clust_2 = c(svg2, rep(NA, max_ln - length(svg2))),
    GCL = c(svgGCL, rep(NA, max_ln - length(svgGCL))),
    SGZ = c(svgSGZ, rep(NA, max_ln - length(svgSGZ))),
    CA4 = c(svgCA4, rep(NA, max_ln - length(svgCA4))),
    CA3 = c(svgCA3, rep(NA, max_ln - length(svgCA3))),
    ML = c(svgML, rep(NA, max_ln - length(svgML))),
    clust_8 = c(svg8, rep(NA, max_ln - length(svg8)))
)

# directory to save .csv list
dir_outputs <- here("processed-data", "nnSVG", "BayesSpace")
fn_out <- file.path(dir_outputs, "DG_BayesSpace_binomial_enriched_SVG_overlapping_markers")

# Export summary as .csv file
write.csv(gene_markers,fn_out, row.names = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
