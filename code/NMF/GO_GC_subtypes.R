###############################
# HPC project
# GSEA of GC.3 & GC.4
# Anthony Ramnauth, Jan 26 2024
###############################

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(sessioninfo)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(rrvgo)
    library(pheatmap)
    library(sessioninfo)
})

load(here::here('sce_enrichment_stats.rda'))

GC.3 <- data.frame(
    ensembl = stats_enrichment$ensembl,
    gene = stats_enrichment$gene,
    logFC = stats_enrichment$logFC_GC.3,
    FDR = stats_enrichment$fdr_GC.3
)

GC.3_entrez <- bitr(GC.3$ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
GC.3_entrez <- GC.3_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

GC.3 <- GC.3[GC.3$ensembl %in% GC.3_entrez$ENSEMBL,]

stopifnot(GC.3$ensembl == GC.3_entrez$ENSEMBL)

GC.3$ENTREZID <- GC.3_entrez$ENTREZID

#filter out FDR > 0.05
GC.3 <- GC.3 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(desc(logFC))

GC.3$group <- "upregulated"
GC.3$group[GC.3$logFC < 0] <- "downregulated"
GC.3$subtype <- "GC3"

GC.4 <- data.frame(
    ensembl = stats_enrichment$ensembl,
    gene = stats_enrichment$gene,
    logFC = stats_enrichment$logFC_GC.4,
    FDR = stats_enrichment$fdr_GC.4
)

GC.4_entrez <- bitr(GC.4$ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 1:many mapping so remove duplicated ENSEMBL mappings
GC.4_entrez <- GC.4_entrez %>%
    distinct(ENSEMBL, .keep_all = TRUE)

GC.4 <- GC.4[GC.4$ensembl %in% GC.4_entrez$ENSEMBL,]

stopifnot(GC.4$ensembl == GC.4_entrez$ENSEMBL)

GC.4$ENTREZID <- GC.4_entrez$ENTREZID

#filter out FDR > 0.05
GC.4 <- GC.4 %>%
    filter(FDR < 0.05) %>%
    dplyr::arrange(desc(logFC))

GC.4$group <- "upregulated"
GC.4$group[GC.4$logFC < 0] <- "downregulated"
GC.4$subtype <- "GC4"

# make the total data frame used for cluster compare
clust_compare <- rbind(GC.3, GC.4)

# Run compareCluster

GO_ALL <- compareCluster(ENTREZID | logFC ~ group+subtype,
    data = clust_compare,
    OrgDb = org.Hs.eg.db,
    fun = gseGO,
    ont = "ALL",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-100,
    pvalueCutoff = 0.001,
    pAdjustMethod = "BH"
)

GO_ALL <- setReadable(GO_ALL, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

save(GO_ALL,
    file = here::here("processed-data", "GC3_GC4_GSEA_new.Rdata")
)

# Make some plots for significant GO terms

# After reviewing, select GO terms that will be use for GO_ALL plot

select_terms <- c(
"neurogenesis",
"neuron differentiation",
"generation of neurons",
"regulation of developmental process",
"cell-cell signaling",
"polysome",
"tissue development",
"cell differentiation",
"collagen-containing extracellular matrix",
"animal organ development",
"synaptic signaling",
"nuclear-transcribed mRNA catabolic process",
"vesicle tethering complex",
"transmembrane transport",
"neuron projection",
"glycosaminoglycan binding",
"spliceosomal complex",
"G protein-coupled receptor signaling pathway",
"positive regulation of developmental process",
"cellular developmental process"
)

pdf(file = here::here("plots", "GC3_GC4_GO_ALL_select.pdf"), width = 5, height = 7)

dotplot(GO_ALL, x="group", showCategory = select_terms) +
    theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(angle = 65, size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
    scale_fill_gradient(low = "black", high = "grey") +
    facet_grid(~factor(subtype, levels=c('GC3', 'GC4')))

dev.off()

pdf(file = here::here("plots", "GC3_GC4_GO_ALL.pdf"), width = 4.5, height = 50)

dotplot(GO_ALL, x="group", showCategory = 1000) +
    theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(angle = 70, size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
    scale_fill_gradient(low = "black", high = "grey") +
    facet_grid(~factor(subtype, levels=c('GC3', 'GC4')))

dev.off()

# cnet plot for genes up and down in GO terms

list <- c("neurogenesis", "tissue development", "neuron differentiation", "animal organ development",
    "cell differentiation", "generation of neurons")

pdf(file = here::here("plots", "GC3_GC4_GO_ALL_cnet.pdf"), width = 80, height = 80)

cnetplot(GO_ALL, list)

dev.off()

pdf(file = here::here("plots", "GC3_GC4_GO_ALL_cnet_catlabel.pdf"))

cnetplot(GO_ALL, list, node_label="category",
        cex_label_category = 1.2)

dev.off()

# directory to save lists
dir_outputs <- here("processed-data")
fn_out_1 <- file.path(dir_outputs, "GO_ALL")

# Export summary as .csv file
write.csv(GO_ALL,fn_out_1, row.names = FALSE)

list_ecm <- c("collagen-containing extracellular matrix", "glycosaminoglycan binding", "collagen trimer")

pdf(file = here::here("plots", "GC3_GC4_GO_ALL_cnet_ECM.pdf"), width = 10, height = 10)

cnetplot(GO_ALL, list_ecm)

dev.off()
