################################################
# spatial_DG_lifespan project
# GO enrichment & plotting from modeling results
# Anthony Ramnauth, May 25 2022
################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(dplyr)
    library(sessioninfo)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(ggplot2)
})

modeling_results <- readRDS(file = here::here("processed-data","pseudobulk_spe","modeling_results.rds"))

enriched_model <- modeling_results$enrichment

## Create named vector of gene lists with ENTREZID (works best) for GO (ranking by t-statistic)

list_1 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_1, t_stat_1)

list_2 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_2, t_stat_2)

list_3 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_3, t_stat_3)

list_4 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_4, t_stat_4)

list_5 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_5, t_stat_5)

list_6 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_6, t_stat_6)

list_7 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_7, t_stat_7)

list_8 <- enriched_model %>%
    dplyr::select(gene, ensembl, fdr_8, t_stat_8)

top_list_1 <- list_1 %>%
    filter(fdr_1 < 0.05) %>%
    dplyr::arrange(desc(t_stat_1)) %>%
    dplyr::select(gene, t_stat_1)

## feature 1: numeric vector
clust_1 = top_list_1[,2]
## feature 2: named vector
names(clust_1) = as.character(top_list_1[,1])
# final vector
clust_1 <- names(clust_1)
clust_1 = bitr(clust_1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   23.08% of input gene IDs are fail to map...

top_list_2 <- list_2 %>%
    filter(fdr_2 < 0.05) %>%
    dplyr::arrange(desc(t_stat_2)) %>%
    dplyr::select(gene, t_stat_2)

## feature 1: numeric vector
clust_2 = top_list_2[,2]
## feature 2: named vector
names(clust_2) = as.character(top_list_2[,1])
# final vector
clust_2 <- names(clust_2)
clust_2 = bitr(clust_2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.46% of input gene IDs are fail to map...

top_list_3 <- list_3 %>%
    filter(fdr_3 < 0.05) %>%
    dplyr::arrange(desc(t_stat_3)) %>%
    dplyr::select(gene, t_stat_3)

## feature 1: numeric vector
GCL = top_list_3[,2]
## feature 2: named vector
names(GCL) = as.character(top_list_3[,1])
# final vector
GCL <- names(GCL)
GCL = bitr(GCL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(GCL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.86% of input gene IDs are fail to map...

top_list_4 <- list_4 %>%
    filter(fdr_4 < 0.05) %>%
    dplyr::arrange(desc(t_stat_4)) %>%
    dplyr::select(gene, t_stat_4)

## feature 1: numeric vector
SGZ = top_list_4[,2]
## feature 2: named vector
names(SGZ) = as.character(top_list_4[,1])
# final vector
SGZ <- names(SGZ)
SGZ = bitr(SGZ, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

top_list_5 <- list_5 %>%
    filter(fdr_5 < 0.05) %>%
    dplyr::arrange(desc(t_stat_5)) %>%
    dplyr::select(gene, t_stat_5)

## feature 1: numeric vector
CA4 = top_list_5[,2]
## feature 2: named vector
names(CA4) = as.character(top_list_5[,1])
# final vector
CA4 <- names(CA4)
CA4 = bitr(CA4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(CA4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   1.79% of input gene IDs are fail to map...

top_list_6 <- list_6 %>%
    filter(fdr_6 < 0.05) %>%
    dplyr::arrange(desc(t_stat_6)) %>%
    dplyr::select(gene, t_stat_6)

## feature 1: numeric vector
CA3 = top_list_6[,2]
## feature 2: named vector
names(CA3) = as.character(top_list_6[,1])
# final vector
CA3 <- names(CA3)
CA3 = bitr(CA3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(CA3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.84% of input gene IDs are fail to map...

list_7 <- list_7 %>%
    arrange(desc(t_stat_7))
# BayesSpace cluster 7 (Molec. Layer) has only 1 gene that is FDR<0.05, so do top 20 instead
top_list_7 <- head(list_7, 20) %>%
    dplyr::select(gene, t_stat_7)

## feature 1: numeric vector
ML = top_list_7[,2]
## feature 2: named vector
names(ML) = as.character(top_list_7[,1])
# final vector
ML <- names(ML)
ML = bitr(ML, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(ML, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   5% of input gene IDs are fail to map...

top_list_8 <- list_8 %>%
    filter(fdr_8 < 0.05) %>%
    dplyr::arrange(desc(t_stat_8)) %>%
    dplyr::select(gene, t_stat_8)

## feature 1: numeric vector
clust_8 = top_list_8[,2]
## feature 2: named vector
names(clust_8) = as.character(top_list_8[,1])
# final vector
clust_8 <- names(clust_8)
clust_8 = bitr(clust_8, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_8, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.48% of input gene IDs are fail to map...

## GO Gene Over Representation Analysis for enrichment in sub-ontologies

clust_1_CC <- enrichGO(gene = clust_1$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_1_MF <- enrichGO(gene = clust_1$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_1_BP <- enrichGO(gene = clust_1$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_2_CC <- enrichGO(gene = clust_2$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_2_MF <- enrichGO(gene = clust_2$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_2_BP <- enrichGO(gene = clust_2$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

GCL_CC <- enrichGO(gene = GCL$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

GCL_MF <- enrichGO(gene = GCL$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

GCL_BP <- enrichGO(gene = GCL$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

SGZ_CC <- enrichGO(gene = SGZ$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

SGZ_MF <- enrichGO(gene = SGZ$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

SGZ_BP <- enrichGO(gene = SGZ$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

CA4_CC <- enrichGO(gene = CA4$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

CA4_MF <- enrichGO(gene = CA4$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

CA4_BP <- enrichGO(gene = CA4$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

CA3_CC <- enrichGO(gene = CA3$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

CA3_MF <- enrichGO(gene = CA3$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

CA3_BP <- enrichGO(gene = CA3$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

ML_CC <- enrichGO(gene = ML$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

ML_MF <- enrichGO(gene = ML$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

ML_BP <- enrichGO(gene = ML$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_8_CC <- enrichGO(gene = clust_8$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_8_MF <- enrichGO(gene = clust_8$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

clust_8_BP <- enrichGO(gene = clust_8$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    readable = TRUE,
    )

## Barplots of subontologies

# For Cluster 1
pdf(file = here::here("plots","pseudobulked","clust_1_GO.pdf"), width = 10, height = 10)

mutate(clust_1_CC, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("cluster_1 Cellular Compartment")
mutate(clust_1_MF, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("cluster_1 Molecular Function")
mutate(clust_1_BP, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("cluster_1 Biological Process")

clust_1_CCx <- setReadable(clust_1_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_1_CCx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("cluster_1 Cellular Compartment Network") +
    theme(legend.position="none")
clust_1_MFx <- setReadable(clust_1_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_1_MFx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("cluster_1 Molecular Function Network") +
    theme(legend.position="none")
clust_1_BPx <- setReadable(clust_1_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_1_BPx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("cluster_1 Biological Process Network") +
    theme(legend.position="none")

clust_1_CCp <- pairwise_termsim(clust_1_CC)
emapplot(clust_1_CCp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("cluster_1 Cellular Compartment Modules")
clust_1_MFp <- pairwise_termsim(clust_1_MF)
emapplot(clust_1_MFp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("cluster_1 Molecular Function Modules")
clust_1_BPp <- pairwise_termsim(clust_1_BP)
emapplot(clust_1_BPp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("cluster_1 Biological Process Modules")

dev.off()

# For Cluster 2
pdf(file = here::here("plots","pseudobulked","clust_2_GO.pdf"), width = 10, height = 10)

mutate(clust_2_CC, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("cluster_2 Cellular Compartment")
mutate(clust_2_MF, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("cluster_2 Molecular Function")
mutate(clust_2_BP, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("cluster_2 Biological Process")

clust_2_CCx <- setReadable(clust_2_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_2_CCx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("cluster_2 Cellular Compartment Network") +
    theme(legend.position="none")
clust_2_MFx <- setReadable(clust_2_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_2_MFx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("cluster_2 Molecular Function Network") +
    theme(legend.position="none")
clust_2_BPx <- setReadable(clust_2_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_2_BPx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("cluster_2 Biological Process Network") +
    theme(legend.position="none")

clust_2_CCp <- pairwise_termsim(clust_2_CC)
emapplot(clust_2_CCp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("cluster_2 Cellular Compartment Modules")
clust_2_MFp <- pairwise_termsim(clust_2_MF)
emapplot(clust_2_MFp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("cluster_2 Molecular Function Modules")
clust_2_BPp <- pairwise_termsim(clust_2_BP)
emapplot(clust_2_BPp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("cluster_2 Biological Process Modules")

dev.off()

# For GCL

pdf(file = here::here("plots","pseudobulked","GCL_GO.pdf"), width = 10, height = 10)

mutate(GCL_CC, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("GCL Cellular Compartment")
mutate(GCL_MF, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("GCL Molecular Function")
mutate(GCL_BP, qscore = -log(p.adjust, base=10)) %>%
    barplot(x="qscore", showCategory=20) +
    ggtitle("GCL Biological Process")

GCL_CCx <- setReadable(GCL_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GCL_CCx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("GCL Cellular Compartment Network") +
    theme(legend.position="none")
GCL_MFx <- setReadable(GCL_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GCL_MFx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("GCL Molecular Function Network") +
    theme(legend.position="none")
GCL_BPx <- setReadable(GCL_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GCL_BPx, showCategory = 20, colorEdge = TRUE, cex_label_category = 0.5,
    cex_label_gene = 0.5) +
    ggtitle("GCL Biological Process Network") +
    theme(legend.position="none")

GCL_CCp <- pairwise_termsim(GCL_CC)
emapplot(GCL_CCp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("GCL Cellular Compartment Modules")
GCL_MFp <- pairwise_termsim(GCL_MF)
emapplot(GCL_MFp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("GCL Molecular Function Modules")
GCL_BPp <- pairwise_termsim(GCL_BP)
emapplot(GCL_BPp, showCategory = 20, color = "p.adjust", cex_label_category = 0.5) +
    ggtitle("GCL Biological Process Modules")

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
