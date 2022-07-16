###############################
# spatial_DG_lifespan project
# Parse DE analysis results
# Anthony Ramnauth, May 12 2022
###############################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(vioplot)
    library(scater)
    library(sessioninfo)
})

# Load modeling results
load(file = here::here("processed-data", "pseudobulk_spe", "infant_bayesspace_cluster_modeling_results.Rdata"))

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

## subset spe data based on age bin
spe_pseudo <- spe_pseudo[, spe_pseudo$age_bin %in% c("Infant")]

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) <- rownames(eb_contrasts)
fdrs0_contrasts <- apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) <- rownames(eb_contrasts)
summary(fdrs0_contrasts < 0.05)

# Merge statistics
f_merge <- function(p, fdr, t) {
    colnames(p) <- paste0("p_value_", colnames(p))
    colnames(fdr) <- paste0("fdr_", colnames(fdr))
    colnames(t) <- paste0("t_stat_", colnames(t))
    res <- as.data.frame(cbind(t, p, fdr))
    res$ensembl <- rownames(res)
    ## Check it's all in order
    stopifnot(identical(rownames(res), rownames(spe_pseudo)))
    res$gene <- rowData(spe_pseudo)$gene_name
    rownames(res) <- NULL
    return(res)
}

results_specificity <-
    f_merge(p = pvals0_contrasts, fdr = fdrs0_contrasts, t = t0_contrasts)
options(width = 400)
head(results_specificity)

pvals_contrasts <- eb_contrasts$p.value
fdrs_contrasts <- apply(pvals_contrasts, 2, p.adjust, "fdr")
dim(pvals_contrasts)
summary(fdrs_contrasts < 0.05)

results_pairwise <-
    f_merge(p = pvals_contrasts, fdr = fdrs_contrasts, t = eb_contrasts$t)
colnames(results_pairwise)
sort(colSums(fdrs_contrasts < 0.05))

f_sig <- function(type, cut = 0.05) {
    cbind(
        "n" = addmargins(table(f_stats[[type]] < cut)),
        "ratio" = addmargins(table(f_stats[[type]] < cut)) / nrow(f_stats)
    )
}
f_sig("full_fdr")

# Match the colnames to the new style
f_rename <- function(x, old, new = old) {
    old_patt <- paste0("_", old, "$")
    i <- grep(old_patt, colnames(x))
    tmp <- gsub(old_patt, "", colnames(x)[i])
    tmp <- paste0(new, "_", tmp)
    colnames(x)[i] <- tmp
    return(x)
}

results_anova <-
    f_rename(f_rename(f_rename(
        f_rename(f_stats, "f", "f_stat"), "p_value"
    ), "fdr"), "Amean")
head(results_anova)

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe", "infant_modeling_results.rds"))

length(which(modeling_results$enrichment$fdr_1 < 0.05))
# 0
length(which(modeling_results$enrichment$fdr_2 < 0.05))
# 0
length(which(modeling_results$enrichment$fdr_3 < 0.05))
# 3
length(which(modeling_results$enrichment$fdr_4 < 0.05))
# 0
length(which(modeling_results$enrichment$fdr_5 < 0.05))
# 0
length(which(modeling_results$enrichment$fdr_6 < 0.05))
# 0
length(which(modeling_results$enrichment$fdr_7 < 0.05))
# 0
length(which(modeling_results$enrichment$fdr_8 < 0.05))
# 2239
## Weird that so many are zero!

cluster <- c(1:8)
genes <- c(0, 0, 3, 0, 0, 0, 0, 2239)
df <- data.frame(cluster, genes)
pdf(file = here::here("plots", "pseudobulked", "plot_infant_enrichment_DEGs.pdf"))
plot(df$genes ~ df$cluster)
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
