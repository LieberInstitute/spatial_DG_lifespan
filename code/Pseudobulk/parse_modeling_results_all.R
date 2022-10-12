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
load(file = here::here("processed-data", "pseudobulk_spe", "bayesspace_cluster_modeling_results.Rdata"))

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

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
#     1               2               3               4               5               6               7
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:11285     FALSE:8021      FALSE:11218     FALSE:7671      FALSE:5735      FALSE:7677      FALSE:10749
# TRUE :8         TRUE :3272      TRUE :75        TRUE :3622      TRUE :5558      TRUE :3616      TRUE :544
#     8
# Mode :logical
# FALSE:11261
# TRUE :32

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
# [1] 14409    66
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
#  f_stat_full p_value_full     fdr_full full_AveExpr         ensembl       gene
#1    19.08982 6.517947e-22 8.839257e-22     2.579431 ENSG00000237491  LINC01409
#2    34.93502 6.524805e-33 1.531453e-32     4.084420 ENSG00000228794  LINC01128
#3    47.26826 3.483369e-39 1.245147e-38     5.043462 ENSG00000188976      NOC2L
#4    20.79932 2.471934e-23 3.561809e-23     2.776608 ENSG00000187961     KLHL17
#5    11.27791 2.593258e-14 2.771361e-14     1.735252 ENSG00000272512 AL645608.7
#6    44.25629 8.827304e-38 2.862120e-37     5.494125 ENSG00000188290       HES4

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

length(which(modeling_results$enrichment$fdr_1 < 0.05))
# [1] 8
length(which(modeling_results$enrichment$fdr_2 < 0.05))
# [1] 3272
length(which(modeling_results$enrichment$fdr_3 < 0.05))
# [1] 75
length(which(modeling_results$enrichment$fdr_4 < 0.05))
# [1] 3622
length(which(modeling_results$enrichment$fdr_5 < 0.05))
# [1] 5558
length(which(modeling_results$enrichment$fdr_6 < 0.05))
# [1] 3616
length(which(modeling_results$enrichment$fdr_7 < 0.05))
# [1] 544
length(which(modeling_results$enrichment$fdr_8 < 0.05))
# [1] 32


cluster <- c(1, 2, 3, 4, 5, 6, 7, 8)
genes <- c(8, 3272, 75, 3622, 5558, 3616, 544, 32)
df <- data.frame(cluster, genes)
pdf(file = here::here("plots", "pseudobulked", "plot_enrichment_DEGs.pdf"))
plot(df$genes ~ df$cluster)
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
