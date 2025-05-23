##################################################
# spatial_DG_lifespan project
# Parse DE analysis results for Manual Annotations
# Anthony Ramnauth, Oct 11 2022
##################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(vioplot)
    library(scater)
    library(ggplot2)
})

# Load modeling results
load(file = here::here("processed-data", "pseudobulk_spe", "manual_annotation_modeling_results.Rdata"))

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "manual_annotated_pseudobulk_spe.rds"))

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
#    CA4              CP             GCL              ML           PCL_CA1         PCL_CA3
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:16835     FALSE:9185      FALSE:9916      FALSE:16659     FALSE:16517     FALSE:14373
# TRUE :144       TRUE :7794      TRUE :7063      TRUE :320       TRUE :462       TRUE :2606
#    SGZ              SL             SLM              SO              SR             SUB
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:16927     FALSE:15782     FALSE:14707     FALSE:15795     FALSE:16751     FALSE:11212
# TRUE :52        TRUE :1197      TRUE :2272      TRUE :1184      TRUE :228       TRUE :5767
#    THAL             WM
# Mode :logical   Mode :logical
# FALSE:16918     FALSE:14057
# TRUE :61        TRUE :2922

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
# [1] 16979   120

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
#1    3.863891 8.511915e-06 9.565412e-06    0.7662687 ENSG00000241860 AL627309.5
#2   13.580893 5.814901e-21 8.156234e-21    2.1969784 ENSG00000237491  LINC01409
#3   84.922590 2.693034e-64 6.407654e-64    4.0752877 ENSG00000228794  LINC01128
#4    1.374014 1.673138e-01 1.694394e-01   -0.2647496 ENSG00000225880  LINC00115
#5    1.038317 4.197731e-01 4.219600e-01    0.2218960 ENSG00000230368     FAM41C
#6    8.657303 5.152916e-14 6.658399e-14    1.0334963 ENSG00000187634     SAMD11

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe",
    "manual_annotated_modeling_results.rds"))
