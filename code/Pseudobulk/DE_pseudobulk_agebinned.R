###################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked BayesSpace clusters
# Anthony Ramnauth, Sept 1 2022
###################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(sessioninfo)
    library(SingleCellExperiment)
    library(rafalib)
    library(limma)
})

# Load SPE
sce_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "First_two_pseudobulk_spe.rds"))

# Format spe object for DE models
colData(sce_pseudo) <- colData(sce_pseudo)[, sort(c(
    "age_bin",
    "BayesSpace",
    "sex",
    "race",
    "ncells"
))]

colData(sce_pseudo)

colData(sce_pseudo)$ncells <- as.numeric(colData(sce_pseudo)$ncells)
colData(sce_pseudo)$race <- as.factor(colData(sce_pseudo)$race)
colData(sce_pseudo)$sex <- as.factor(colData(sce_pseudo)$sex)

colData(sce_pseudo)

# Extract the data
mat <- assays(sce_pseudo)$logcounts

# Matrix for a regression-like model
mod <- with(colData(sce_pseudo), model.matrix(~ 0 + BayesSpace))

# Compute correlation
corfit <- duplicateCorrelation(mat, mod, block = sce_pseudo$age_bin)

######### ENRICHMENT t-stats ######################

cluster_id <- splitit(sce_pseudo$BayesSpace)

eb0_list <- lapply(cluster_id, function(x) {
    res <- rep(0, ncol(sce_pseudo))
    res[x] <- 1
    m <- model.matrix(~res)
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudo$age_bin,
            correlation = corfit$consensus.correlation
        )
    )
})

######### PAIRWISE t-stats ######################

fit <-
    lmFit(
        mat,
        design = mod,
        block = sce_pseudo$age_bin,
        correlation = corfit$consensus.correlation
    )
eb <- eBayes(fit)

## Define the contrasts for each group vs another one
cluster_combs <- combn(colnames(mod), 2)
cluster_constrats <- apply(cluster_combs, 2, function(x) {
    z <- paste(x, collapse = "-")
    makeContrasts(contrasts = z, levels = mod)
})
rownames(cluster_constrats) <- colnames(mod)
colnames(cluster_constrats) <-
    apply(cluster_combs, 2, paste, collapse = "-")
eb_contrasts <- eBayes(contrasts.fit(fit, cluster_constrats))

######### ANOVA t-stats ######################

## From layer_specificity.R
fit_f_model <- function(sce_pseudo) {

    ## Extract the data
    mat <- assays(sce_pseudo)$logcounts


    ## Build a group model
    mod <- with(colData(sce_pseudo), model.matrix(~BayesSpace))
    colnames(mod) <- gsub("clusters", "", colnames(mod))

    # already made in beginning of script #remember to adjust for age or sex
    corfit <-
        duplicateCorrelation(mat, mod, block = sce_pseudo$age_bin)
    message(paste(Sys.time(), "correlation:", corfit$consensus.correlation))
    fit <-
        lmFit(
            mat,
            design = mod,
            block = sce_pseudo$age_bin,
            correlation = corfit$consensus.correlation
        )
    eb <- eBayes(fit)
    return(eb)
}

ebF_list <-
    lapply(list("full" = sce_pseudo), fit_f_model)

## Extract F-statistics
f_stats <- do.call(cbind, lapply(names(ebF_list), function(i) {
    x <- ebF_list[[i]]
    top <-
        topTable(
            x,
            coef = 2:8,
            sort.by = "none",
            number = length(x$F)
        )
    # identical(p.adjust(top$P.Value, 'fdr'), top$adj.P.Val)
    res <- data.frame(
        "f" = top$F,
        "p_value" = top$P.Value,
        "fdr" = top$adj.P.Val,
        "AveExpr" = top$AveExpr,
        stringsAsFactors = FALSE
    )
    colnames(res) <- paste0(i, "_", colnames(res))
    return(res)
}))
f_stats$ensembl <- rownames(sce_pseudo)
f_stats$gene <- rowData(sce_pseudo)$gene_name
rownames(f_stats) <- NULL

head(f_stats)

save(f_stats, eb0_list, eb_contrasts,
    file = here::here("processed-data", "pseudobulk_spe", "slide1_2_bayesspace_cluster_modeling_results.Rdata")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
