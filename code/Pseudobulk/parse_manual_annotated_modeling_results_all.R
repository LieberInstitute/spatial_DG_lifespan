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
    library(sessioninfo)
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
#    CA4              CP             GCL              ML           PCL_CA1         PCL_CA3           SGZ
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:14557     FALSE:8925      FALSE:10708     FALSE:14591     FALSE:14468     FALSE:13708     FALSE:14591
# TRUE :34        TRUE :5666      TRUE :3883                      TRUE :123       TRUE :883
#     SL             SLM              SO              SR             SUB              WM
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:14393     FALSE:14580     FALSE:14063     FALSE:14583     FALSE:14591     FALSE:12916
# TRUE :198       TRUE :11        TRUE :528       TRUE :8                         TRUE :1675

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
#  f_stat_full p_value_full     fdr_full full_AveExpr         ensembl      gene
#1   11.096881 5.450360e-13 6.488228e-13     2.346724 ENSG00000237491 LINC01409
#2   66.603684 1.332830e-35 2.925289e-35     4.085475 ENSG00000228794 LINC01128
#3    3.912202 5.989356e-05 6.186952e-05     0.717042 ENSG00000225880 LINC00115
#4   11.372171 3.100797e-13 3.713373e-13     1.665840 ENSG00000187634    SAMD11
#5  239.134296 9.000950e-55 4.260971e-54     5.338135 ENSG00000188976     NOC2L
#6   21.421764 2.587459e-20 3.641359e-20     2.733696 ENSG00000187961    KLHL17

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "pseudobulk_spe",
    "manual_annotated_modeling_results.rds"))

length(which(modeling_results$enrichment$fdr_CA4 < 0.05))
# [1] 34
length(which(modeling_results$enrichment$fdr_CP < 0.05))
# [1] 5666
length(which(modeling_results$enrichment$fdr_GCL < 0.05))
# [1] 3883
length(which(modeling_results$enrichment$fdr_ML < 0.05))
# [1] 0
length(which(modeling_results$enrichment$fdr_PCL_CA1 < 0.05))
# [1] 123
length(which(modeling_results$enrichment$fdr_PCL_CA3 < 0.05))
# [1] 883
length(which(modeling_results$enrichment$fdr_SGZ < 0.05))
# [1] 0
length(which(modeling_results$enrichment$fdr_SL < 0.05))
# [1] 198
length(which(modeling_results$enrichment$fdr_SLM < 0.05))
# [1] 11
length(which(modeling_results$enrichment$fdr_SO < 0.05))
# [1] 528
length(which(modeling_results$enrichment$fdr_SR < 0.05))
# [1] 8
length(which(modeling_results$enrichment$fdr_SUB < 0.05))
# [1] 0
length(which(modeling_results$enrichment$fdr_WM < 0.05))
# [1] 1675

cluster <- c("CA4" = 1, "CP" = 2, "GCL" = 3, "ML" = 4, "PCL_CA1" = 5,
    "PCL_CA3" = 6, "SGZ" = 7, "SL" = 8, "SLM" = 9, "SO" = 10, "SR" = 11, "SUB" = 12, "WM" = 13)
genes <- c(34, 5666, 3883, 0, 123, 883, 0, 198, 11, 528, 8, 0, 1675)
df <- data.frame(cluster, genes)
pdf(file = here::here("plots", "pseudobulked", "plot_manual_annotated_enrichment_DEGs.pdf"))
plot(df$genes ~ df$cluster, xaxt = "n")
axis(1,
    at = 1:13,
     labels = c("CA4", "CP", "GCL", "ML", "CA1", "CA3", "SGZ", "SL", "SLM",
         "SO", "SR", "SUB", "WM"))
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
