###############################
# spatial_DG_lifespan project
# Parse DE analysis results
# Anthony Ramnauth, May 12 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
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
#     1               2               3               4
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:10532     FALSE:12864     FALSE:5019      FALSE:7961
# TRUE :2823      TRUE :491       TRUE :8336      TRUE :5394
#     5               6               7               8
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:12640     FALSE:12922     FALSE:6279      FALSE:12979
# TRUE :715       TRUE :433       TRUE :7076      TRUE :376
#     9               10
# Mode :logical   Mode :logical
# FALSE:9825      FALSE:7708
# TRUE :3530      TRUE :5647

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
# [1] 13355    66
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
#          n        ratio
#FALSE     7 0.0005241483
#TRUE  13348 0.9994758517
#Sum   13355 1.0000000000

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
#  f_stat_full  p_value_full      fdr_full full_AveExpr         ensembl      gene
#1   19.660516  1.288494e-23  1.476181e-23     2.278726 ENSG00000237491 LINC01409
#2  122.918068  3.779667e-67  6.958567e-67     4.021929 ENSG00000228794 LINC01128
#3    8.775102  8.781450e-12  9.045604e-12     1.190339 ENSG00000187634    SAMD11
#4  693.711329 3.107382e-118 1.361072e-117     5.239307 ENSG00000188976     NOC2L
#5   20.666296  1.603342e-24  1.851343e-24     2.411745 ENSG00000187961    KLHL17
#6  366.775717  5.827790e-99  1.781418e-98     5.686110 ENSG00000188290      HES4

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

length(which(modeling_results$enrichment$fdr_1 < 0.05))
# [1] 2823
length(which(modeling_results$enrichment$fdr_2 < 0.05))
# [1] 491
length(which(modeling_results$enrichment$fdr_3 < 0.05))
# [1] 8336
length(which(modeling_results$enrichment$fdr_4 < 0.05))
# [1] 5394
length(which(modeling_results$enrichment$fdr_6 < 0.05))
# [1] 433
length(which(modeling_results$enrichment$fdr_7 < 0.05))
# [1] 7076
length(which(modeling_results$enrichment$fdr_8 < 0.05))
# [1] 376
length(which(modeling_results$enrichment$fdr_9 < 0.05))
# [1] 3530
length(which(modeling_results$enrichment$fdr_10 < 0.05))
# [1] 5647

cluster <- c(1, 2, 3, 4, 6, 7, 8, 9, 10)
genes <- c(2823, 491, 8336, 5394, 433, 7076, 376, 3530, 5647)
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
