###################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked BayesSpace clusters
# Anthony Ramnauth, May 27 2022
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
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on age bin
noinfant_spe_pseudo <- spe_pseudo[, !spe_pseudo$age_bin %in% c("Infant")]

# Format spe object for DE models
colData(noinfant_spe_pseudo) <- colData(noinfant_spe_pseudo)[, sort(c(
    "age",
    "sample_id",
    "BayesSpace",
    "sex",
    "race",
    "ncells"
))]

colData(noinfant_spe_pseudo)

colData(noinfant_spe_pseudo)$ncells <- as.numeric(colData(noinfant_spe_pseudo)$ncells)
colData(noinfant_spe_pseudo)$race <- as.factor(colData(noinfant_spe_pseudo)$race)
colData(noinfant_spe_pseudo)$sample_id <- as.factor(colData(noinfant_spe_pseudo)$sample_id)
colData(noinfant_spe_pseudo)$sex <- as.factor(colData(noinfant_spe_pseudo)$sex)

colData(noinfant_spe_pseudo)

# Drop things we don't need
spatialCoords(noinfant_spe_pseudo) <- NULL
imgData(noinfant_spe_pseudo) <- NULL

# Extract the data
mat <- assays(noinfant_spe_pseudo)$logcounts

# Make mat_formula
var_i <- "BayesSpace"
covars <- c("age", "sex")
mat_formula <- eval(str2expression(paste("~", "0", "+", var_i, "+", paste(covars, collapse = " + "))))

# Matrix for a regression-like model
mod <- model.matrix(mat_formula, data = colData(noinfant_spe_pseudo))

# Compute correlation
corfit <- duplicateCorrelation(mat, mod, block = noinfant_spe_pseudo$sample_id)

######### ENRICHMENT t-stats ######################

cluster_id <- splitit(colData(noinfant_spe_pseudo)[, var_i])

eb0_list <- lapply(cluster_id, function(x) {
    res <- rep(0, ncol(noinfant_spe_pseudo))
    res[x] <- 1
    res_formula <- paste("~", "res", "+", paste(covars, collapse = " + "))
    m <- with(
        colData(noinfant_spe_pseudo),
        model.matrix(eval(str2expression(res_formula)))
    )
    eBayes(
        lmFit(
            mat,
            design = m,
            block = noinfant_spe_pseudo$sample_id,
            correlation = corfit$consensus.correlation
        )
    )
})


######### PAIRWISE t-stats ######################

fit <-
    lmFit(
        mat,
        design = mod,
        block = noinfant_spe_pseudo$sample_id,
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
fit_f_model <- function(spe) {

    ## Extract the data
    mat <- assays(spe)$logcounts

    ## For dropping un-used levels
    colData(spe)[[var_i]] <- as.factor(colData(spe)[[var_i]])

    ## Build a group model
    # already made in beginning of script #remember to adjust for age or sex
    corfit <-
        duplicateCorrelation(mat, mod, block = spe$sample_id)
    message(paste(Sys.time(), "correlation:", corfit$consensus.correlation))
    fit <-
        lmFit(
            mat,
            design = mod,
            block = spe$sample_id,
            correlation = corfit$consensus.correlation
        )
    eb <- eBayes(fit)
    return(eb)
}

ebF_list <-
    lapply(list("full" = noinfant_spe_pseudo), fit_f_model)

## Extract F-statistics
f_stats <- do.call(cbind, lapply(names(ebF_list), function(i) {
    x <- ebF_list[[i]]
    top <-
        topTable(
            x,
            coef = 2:ncol(x$coefficients),
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
f_stats$ensembl <- rownames(noinfant_spe_pseudo)
f_stats$gene <- rowData(noinfant_spe_pseudo)$gene_name
rownames(f_stats) <- NULL

head(f_stats)

save(f_stats, eb0_list, eb_contrasts,
    file = here::here("processed-data", "pseudobulk_spe", "NOinfant_bayesspace_cluster_modeling_results.Rdata")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
