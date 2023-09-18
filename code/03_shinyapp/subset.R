setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

library(spatialLIBD)
library(lobstr)
library(dplyr)
library(here)

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# Order by age
spe <- spe[, order(spe$age)]

## Check how big it is in memory
lobstr::obj_size(spe)
# 4.93 GB GB
## That's too big for shinyapps.io. Aim to have an object near 2GB.

## Subset the spe object outside of shinyapps.io. Otherwise, the peak memory is
## still affected by loading the object.
## Also, running lobstr::obj_size() takes a while to run, which we don't need
## to run every time someone accesses the shiny app.
imgData(spe) <-
    imgData(spe)[!imgData(spe)$image_id %in% c("hires", "detected", "aligned"), ]
assays(spe)$counts <- NULL
lobstr::obj_size(spe)
# 1.79 GB
## Ok, this is reasonable.

# Remove columns in colData(spe) that are unnecessary for Shiny app

# Don't need 10x clusters
colData(spe)[, c(6:15)] <- NULL
# Don't need rin
spe$rin <- NULL
# Don't need NBW to sizeFactor
colData(spe)[, c(17:31)] <- NULL
#Don't need in_tissue
spe$in_tissue <- NULL

# Rename bayesSpace_harmony_10 to BayesSpace
df <-
    data.frame(spe$key, spe$sample_id, spe$bayesSpace_harmony_10)

df <- df %>%
    mutate(
        BayesSpace = case_when(
            spe.bayesSpace_harmony_10 == 1 ~ "SLM",
            spe.bayesSpace_harmony_10 == 2 ~ "ML",
            spe.bayesSpace_harmony_10 == 3 ~ "CP",
            spe.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe.bayesSpace_harmony_10 == 5 ~ "SR",
            spe.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe.bayesSpace_harmony_10 == 7 ~ "GCL",
            spe.bayesSpace_harmony_10 == 8 ~ "SL",
            spe.bayesSpace_harmony_10 == 9 ~ "CA1",
            spe.bayesSpace_harmony_10 == 10 ~ "WM",
        )
    )

colData(spe)$BayesSpace <-
    factor(df$BayesSpace)

spe$bayesSpace_harmony_10 <- NULL

## Save the reduced version of the spe object in the shiny app directory
## instead of using soft links.
saveRDS(spe, file = here::here("code", "03_shinyapp", "spe.rds"))

# Load spe_pseudo object
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

lobstr::obj_size(spe_pseudo)
# 1.60 GB

# Give BayesSpace clusters anatomically recognized names
df1 <-
    data.frame(spe_pseudo$key, spe_pseudo$sample_id, spe_pseudo$BayesSpace)

df1 <- df1 %>%
    mutate(
        BayesSpace = case_when(
            spe_pseudo.BayesSpace == 1 ~ "SLM",
            spe_pseudo.BayesSpace == 2 ~ "ML",
            spe_pseudo.BayesSpace == 4 ~ "CA3_4",
            spe_pseudo.BayesSpace == 5 ~ "SR",
            spe_pseudo.BayesSpace == 6 ~ "SGZ",
            spe_pseudo.BayesSpace == 7 ~ "GCL",
            spe_pseudo.BayesSpace == 8 ~ "SL",
            spe_pseudo.BayesSpace == 9 ~ "CA1",
            spe_pseudo.BayesSpace == 10 ~ "WM",
        )
    )

colData(spe_pseudo)$BayesSpace <- factor(df1$BayesSpace)
spe_pseudo$bayesSpace_harmony_10 <- NULL

## Save the reduced version of the spe object in the shiny app directory
## instead of using soft links.
saveRDS(spe_pseudo, file = here::here("code", "03_shinyapp", "pseudobulk_spe.rds"))

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data", "pseudobulk_spe", "spatialLIBD_modeling_results.rds"))

lobstr::obj_size(modeling_results)
# 21.26 MB

# Save modeling results to Shiny directory
saveRDS(modeling_results, file = here::here("code", "03_shinyapp", "modeling_results.rds"))

## For sig_genes_extract_all() to work https://github.com/LieberInstitute/Visium_IF_AD/blob/5e3518a9d379e90f593f5826cc24ec958f81f4aa/code/05_deploy_app_wholegenome/app.R#L37-L44
spe_pseudo$spatialLIBD <- spe_pseudo$BayesSpace

## Check that we have the right number of tests
k <- 9
tests <- lapply(modeling_results, function(x) {
    colnames(x)[grep("stat", colnames(x))]
})
stopifnot(length(tests$anova) == 1) ## assuming only "all"
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))

sig_genes <- sig_genes_extract_all(
    n = nrow(spe_pseudo),
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

## Check that we have the right number of tests.
## the + 1 at the end assumes only "all"
stopifnot(length(unique(sig_genes$test)) == choose(k, 2) * 2 + k + 1)

lobstr::obj_size(sig_genes)
# 452.61 MB

dim(sig_genes)
# 1051076      13

## Drop parts we don't need to reduce the memory
sig_genes$in_rows <- NULL
sig_genes$in_rows_top20 <- NULL
lobstr::obj_size(sig_genes)
# 91.03 MB

# ## Subset sig_genes
sig_genes <- subset(sig_genes, fdr < 0.05)
dim(sig_genes)
# # [1] 426010     11
lobstr::obj_size(sig_genes)
# # 38.20 MB

# Save subset of sig genes to Shiny directory
saveRDS(sig_genes, file = here::here("code", "03_shinyapp", "sig_genes_subset.rds"))
