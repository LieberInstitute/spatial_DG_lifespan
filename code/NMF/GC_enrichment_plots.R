suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(sessioninfo)
    library(EnhancedVolcano)
    library(dplyr)
    library(RColorBrewer)
})

load(here::here('sce_enrichment_stats.rda'))

GC.3 <- data.frame(
    gene_id = stats_enrichment$ensembl,
    gene_name = stats_enrichment$gene,
    logFC = stats_enrichment$logFC_GC.3,
    FDR = stats_enrichment$fdr_GC.3,
    tstats = stats_enrichment$t_stat_GC.3
)

GC.2 <- data.frame(
    gene_id = stats_enrichment$ensembl,
    gene_name = stats_enrichment$gene,
    logFC = stats_enrichment$logFC_GC.2,
    FDR = stats_enrichment$fdr_GC.2,
    tstats = stats_enrichment$t_stat_GC.2
)

GC.4 <- data.frame(
    gene_id = stats_enrichment$ensembl,
    gene_name = stats_enrichment$gene,
    logFC = stats_enrichment$logFC_GC.4,
    FDR = stats_enrichment$fdr_GC.4,
    tstats = stats_enrichment$t_stat_GC.4
)

## Colors for the significant and not significant genes
keyvals_3 <- ifelse(
    (GC.3$FDR < 0.05) &
        (GC.3$logFC > 1.5 | GC.3$logFC < -1.5), "red", "gray47"
)

## Legend names
names(keyvals_3)[keyvals_3 == "red"] <- "Adjusted P-value < 0.05 & LFC > 1.5"
names(keyvals_3)[keyvals_3 == "gray47"] <- "Not significant"

## Colors for the significant and not significant genes
keyvals_4 <- ifelse(
    (GC.4$FDR < 0.05) &
        (GC.4$logFC > 1.5 | GC.4$logFC < -1.5), "red", "gray47"
)

## Legend names
names(keyvals_4)[keyvals_4 == "red"] <- "Adjusted P-value < 0.05 & LFC > 1.5"
names(keyvals_4)[keyvals_4 == "gray47"] <- "Not significant"

pdf(file = here::here("plots", "GC_enrichment.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(GC.3,
    lab = GC.3$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1.5,
    pCutoff = 0.049,
    max.overlaps = 26,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    colCustom = keyvals_3,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "GC.3 enrichemnt",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-8, 8)) +
    ylim(c(0, 130))

EnhancedVolcano(GC.2,
    lab = GC.2$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1.5,
    pCutoff = 0.049,
    max.overlaps = 10,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "GC.2 enrichemnt"
    )

EnhancedVolcano(GC.4,
    lab = GC.4$gene_name,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1.5,
    pCutoff = 0.049,
    max.overlaps = 26,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    colCustom = keyvals_4,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "GC.4 enrichemnt",
    gridlines.major = FALSE,
    gridlines.minor = FALSE
    ) +
    xlim(c(-8, 8)) +
    ylim(c(0, 130))

dev.off()

# plot tstat vs. tstat for GC.3 and GC.4 since they are on the opposite ends of the spectrum

## colors
cols = rep(1, length(stats_enrichment$enrichment))
cols[stats_enrichment$fdr_GC.3 > 0.05 & stats_enrichment$fdr_GC.4 < 0.05] = 2
cols[stats_enrichment$fdr_GC.3 < 0.05 & stats_enrichment$fdr_GC.4 > 0.05] = 3
cols[stats_enrichment$fdr_GC.3 < 0.05 & stats_enrichment$fdr_GC.4 < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "GC.4"
names(cols)[cols==3] = "GC.3"
names(cols)[cols==4] = "Both"

pdf(file = here::here("plots","GC4_vs_GC3_tstatsneg4.pdf"))

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
	ylim = c(-12, 16), xlim = c(-12, 21), pch = 21,bg=cols,
	xlab = "GC.3 (t-stat)", ylab = "GC.4 (t-stat)",
	main = "T-statistics")
text(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
    labels = ifelse((stats_enrichment$t_stat_GC.4 < -8 & stats_enrichment$fdr_GC.3 > 0.05),
        stats_enrichment$gene, NA), cex = 0.7, pos = 4)
legend("bottomright", c("GC.4", "GC.3", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()

pdf(file = here::here("plots","GC4_vs_GC3_tstatspos4.pdf"))

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
	ylim = c(-12, 16), xlim = c(-12, 21), pch = 21,bg=cols,
	xlab = "GC.3 (t-stat)", ylab = "GC.4 (t-stat)",
	main = "T-statistics")
text(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
    labels = ifelse((stats_enrichment$t_stat_GC.4 > 9.5 & stats_enrichment$fdr_GC.3 > 0.05),
        stats_enrichment$gene, NA), cex = 0.7, pos = 2)
legend("bottomright", c("GC.4", "GC.3", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()

pdf(file = here::here("plots","GC4_vs_GC3_tstatspos3.pdf"))

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
	ylim = c(-12, 16), xlim = c(-12, 21), pch = 21,bg=cols,
	xlab = "GC.3 (t-stat)", ylab = "GC.4 (t-stat)",
	main = "T-statistics")
text(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
    labels = ifelse((stats_enrichment$t_stat_GC.3 > 12 & stats_enrichment$fdr_GC.4 > 0.05),
        stats_enrichment$gene, NA), cex = 0.7, pos = 4)
legend("bottomright", c("GC.4", "GC.3", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()

pdf(file = here::here("plots","GC4_vs_GC3_tstatsneg3.pdf"))

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
	ylim = c(-12, 16), xlim = c(-12, 21), pch = 21,bg=cols,
	xlab = "GC.3 (t-stat)", ylab = "GC.4 (t-stat)",
	main = "T-statistics")
text(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
    labels = ifelse((stats_enrichment$t_stat_GC.3 < -7.5 & stats_enrichment$fdr_GC.4 > 0.05),
        stats_enrichment$gene, NA), cex = 0.7, pos = 3)
legend("bottomright", c("GC.4", "GC.3", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()

pdf(file = here::here("plots","GC4_vs_GC3_tstatsanti4.pdf"))

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
	ylim = c(-12, 16), xlim = c(-12, 21), pch = 21,bg=cols,
	xlab = "GC.3 (t-stat)", ylab = "GC.4 (t-stat)",
	main = "T-statistics")
text(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
    labels = ifelse((stats_enrichment$t_stat_GC.3 > 5 & stats_enrichment$t_stat_GC.4 < -2.5),
        stats_enrichment$gene, NA), cex = 0.7, pos = 4)
legend("bottomright", c("GC.4", "GC.3", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()

pdf(file = here::here("plots","GC4_vs_GC3_tstatsanti3.pdf"))

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
	ylim = c(-12, 16), xlim = c(-12, 21), pch = 21,bg=cols,
	xlab = "GC.3 (t-stat)", ylab = "GC.4 (t-stat)",
	main = "T-statistics")
text(stats_enrichment$t_stat_GC.3, stats_enrichment$t_stat_GC.4,
    labels = ifelse((stats_enrichment$t_stat_GC.4 > 4 & stats_enrichment$t_stat_GC.3 < -3.5),
        stats_enrichment$gene, NA), cex = 0.7, pos = 2)
legend("bottomright", c("GC.4", "GC.3", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()

# Need to derive the numbers of DEGs shared, not shared, etc.

# GC.4 DEGs
table(stats_enrichment$fdr_GC.3 > 0.05 & stats_enrichment$fdr_GC.4 < 0.05)
#FALSE  TRUE
#18627  2279

# GC.3 DEGs
table(stats_enrichment$fdr_GC.3 < 0.05 & stats_enrichment$fdr_GC.4 > 0.05)
#FALSE  TRUE
#18015  2891

#both
table(stats_enrichment$fdr_GC.3 < 0.05 & stats_enrichment$fdr_GC.4 < 0.05)
#FALSE  TRUE
#17957  2949

# get significant DEGs for table
# directory to save  results
dir_outputs <- here("processed-data")

GC.3 <- GC.3 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(FDR)

fn_out1 <- file.path(dir_outputs, "GC3_DE")
# Export summary as .csv file
write.csv(GC.3, fn_out1, row.names = FALSE)

GC.4 <- GC.4 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(FDR)

fn_out <- file.path(dir_outputs, "GC4_DE")
# Export summary as .csv file
write.csv(GC.4, fn_out, row.names = FALSE)


