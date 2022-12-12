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
#     1               2               3               4               6               7               8
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:11528     FALSE:7878      FALSE:10316     FALSE:7266      FALSE:6002      FALSE:10784     FALSE:11531
# TRUE :75        TRUE :3725      TRUE :1287      TRUE :4337      TRUE :5601      TRUE :819       TRUE :72

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
# [1] 11603    36
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
#  f_stat_full p_value_full     fdr_full full_AveExpr         ensembl      gene
#1    14.46728 9.274405e-11 9.400797e-11     2.437106 ENSG00000237491 LINC01409
#2    81.20002 2.697810e-26 3.827202e-26     4.161927 ENSG00000228794 LINC01128
#3   332.83319 3.220475e-41 8.858979e-41     5.359479 ENSG00000188976     NOC2L
#4    20.15705 2.546744e-13 2.643102e-13     2.726867 ENSG00000187961    KLHL17
#5   267.80964 7.221771e-39 1.695553e-38     5.620692 ENSG00000188290      HES4
#6   143.60117 3.267569e-32 5.572252e-32     5.485435 ENSG00000187608     ISG15

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

length(which(modeling_results$enrichment$fdr_1 < 0.05))
# [1] 75
length(which(modeling_results$enrichment$fdr_2 < 0.05))
# [1] 3725
length(which(modeling_results$enrichment$fdr_3 < 0.05))
# [1] 1287
length(which(modeling_results$enrichment$fdr_4 < 0.05))
# [1] 4337
length(which(modeling_results$enrichment$fdr_6 < 0.05))
# [1] 5601
length(which(modeling_results$enrichment$fdr_7 < 0.05))
# [1] 819
length(which(modeling_results$enrichment$fdr_8 < 0.05))
# [1] 72


cluster <- c(1, 2, 3, 4, 6, 7, 8)
genes <- c(75, 3725, 1287, 4337, 5601, 819, 72)
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
