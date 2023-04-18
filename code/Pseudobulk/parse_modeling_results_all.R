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
#     1               2               4               5
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:7875      FALSE:11755     FALSE:6758      FALSE:9961
# TRUE :4943      TRUE :1063      TRUE :6060      TRUE :2857
#     6               7               8               9
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:12116     FALSE:5268      FALSE:11614     FALSE:8454
# TRUE :702       TRUE :7550      TRUE :1204      TRUE :4364
#     10
# Mode :logical
# FALSE:5724
# TRUE :7094

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
# [1] 12818    55
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
#          n       ratio
#FALSE     7 0.000546107
#TRUE  12811 0.999453893
#Sum   12818 1.000000000

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
#1    17.53043  4.798732e-20  5.190292e-20     2.256768 ENSG00000237491 LINC01409
#2   103.17267  6.811564e-58  1.043138e-57     4.058748 ENSG00000228794 LINC01128
#3   667.12879 4.677064e-109 1.615484e-108     5.296257 ENSG00000188976     NOC2L
#4    22.02403  7.776137e-24  8.694567e-24     2.490440 ENSG00000187961    KLHL17
#5   534.19594 9.669470e-103 2.814973e-102     5.682124 ENSG00000188290      HES4
#6   404.25054  7.149063e-95  1.753477e-94     5.432274 ENSG00000187608     ISG15

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe", "modeling_results.rds"))

length(which(modeling_results$enrichment$fdr_1 < 0.05))
# [1] 4943
length(which(modeling_results$enrichment$fdr_2 < 0.05))
# [1] 1063
length(which(modeling_results$enrichment$fdr_4 < 0.05))
# [1] 6060
length(which(modeling_results$enrichment$fdr_6 < 0.05))
# [1] 702
length(which(modeling_results$enrichment$fdr_7 < 0.05))
# [1] 7550
length(which(modeling_results$enrichment$fdr_8 < 0.05))
# [1] 1204
length(which(modeling_results$enrichment$fdr_9 < 0.05))
# [1] 4364
length(which(modeling_results$enrichment$fdr_10 < 0.05))
# [1] 7094

cluster <- c(1, 2, 4, 6, 7, 8, 9, 10)
genes <- c(4943, 1063, 6060, 702, 7550, 1204, 4364, 7094)
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
