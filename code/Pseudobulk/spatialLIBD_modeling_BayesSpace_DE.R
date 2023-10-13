###################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked BayesSpace clusters
# Anthony Ramnauth, Sept 18 2023
###################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(spatialLIBD)
    library(dplyr)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# Give BayesSpace clusters anatomically recognized names
df <-
    data.frame(spe_pseudo$key, spe_pseudo$sample_id, spe_pseudo$BayesSpace)

df <- df %>%
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

colData(spe_pseudo)$BayesSpace <- factor(df$BayesSpace)

## To avoid having to change parameters later on
spe_pseudo$registration_variable <- spe_pseudo$BayesSpace
spe_pseudo$registration_sample_id <- spe_pseudo$sample_id

## Set arguments used in spatialLIBD::registration_wrapper()
covars <- c("age", "sex")
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"

## Taken from spatialLIBD::registration_wrapper()
## https://github.com/LieberInstitute/spatialLIBD/blob/master/R/registration_wrapper.R
registration_mod <-
    registration_model(spe_pseudo, covars = covars)

block_cor <-
    registration_block_cor(spe_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
results_pairwise <-
    registration_stats_pairwise(
        spe_pseudo,
        registration_model = registration_mod,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )

results_anova <-
    registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name,
        suffix = suffix
    )

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

