################################################
# Spatial_DG_lifespan project
# Finding overlapping markers from various tests
# Anthony Ramnauth, July 19 2022
################################################

library(dplyr)
library(sessioninfo)

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
bien1 <- en1[bi1 %in% en1]
svg1 <- svg_list[bien1 %in% svg_list]
svg1 <- svg1[!is.na(svg1)]

# Create vectors of the gene names
en2 <- c(enriched2$gene_name)
bi2 <- c(binom2$gene_name)

# Find gene names from one vector in the other
bien2 <- en2[bi2 %in% en2]
svg2 <- svg_list[bien2 %in% svg_list]
svg2 <- svg2[!is.na(svg2)]

# Create vectors of the gene names
enGCL <- c(enrichedGCL$gene_name)
biGCL <- c(binomGCL$gene_name)

# Find gene names from one vector in the other
bienGCL <- enGCL[biGCL %in% enGCL]
svgGCL <- svg_list[bienGCL %in% svg_list]
svgGCL <- svgGCL[!is.na(svgGCL)]

# Create vectors of the gene names
enSGZ <- c(enrichedSGZ$gene_name)
biSGZ <- c(binomSGZ$gene_name)

# Find gene names from one vector in the other
bienSGZ <- enSGZ[biSGZ %in% enSGZ]
svgSGZ <- svg_list[bienSGZ %in% svg_list]
svgSGZ <- svgSGZ[!is.na(svgSGZ)]

# Create vectors of the gene names
enCA4 <- c(enrichedCA4$gene_name)
biCA4 <- c(binomCA4$gene_name)

# Find gene names from one vector in the other
bienCA4 <- enCA4[biCA4 %in% enCA4]
svgCA4 <- svg_list[bienCA4 %in% svg_list]
svgCA4 <- svgCA4[!is.na(svgCA4)]

# Create vectors of the gene names
enCA3 <- c(enrichedCA3$gene_name)
biCA3 <- c(binomCA3$gene_name)

# Find gene names from one vector in the other
bienCA3 <- enCA3[biCA3 %in% enCA3]
svgCA3 <- svg_list[bienCA3 %in% svg_list]
svgCA3 <- svgCA3[!is.na(svgCA3)]

# Create vectors of the gene names
enML <- c(enrichedML$gene_name)
biML <- c(binomML$gene_name)

# Find gene names from one vector in the other
bienML <- enML[biML %in% enML]
svgML <- svg_list[bienML %in% svg_list]
svgML <- svgML[!is.na(svgML)]

# Create vectors of the gene names
en8 <- c(enriched8$gene_name)
bi8 <- c(binom8$gene_name)

# Find gene names from one vector in the other
bien8 <- en8[bi8 %in% en8]
svg8 <- svg_list[bien8 %in% svg_list]
svg8 <- svg8[!is.na(svg8)]

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
