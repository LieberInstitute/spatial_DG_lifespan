#####################################################
# spatial_HPC project
# NMF pattern projection from human to animal snRNAseq
# Anthony Ramnauth, Dec 07 2023
#####################################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(ggplot2)
    library(scater)
    library(scran)
    library(scry)
    library(CoGAPS)
    library(RcppML)
    library(schex)
    library(reshape2)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(tidyr)
})

# Load mouse DG SCE from Erik Nelson's publication
# https://doi.org/10.1002/hipo.23548
load(file=here::here('sce_subset_mouse.rda'))

set.seed(1000)
sce.subset <- logNormCounts(sce.subset)

# Adatpted from Erik's script for ECS project for PCA and UMAP for plotting purposes
# https://github.com/Erik-D-Nelson/ARG_HPC_snRNAseq/blob/main/code/02_processing/build_sce_QC_processing.R
##re-run analysis from devianceFeatureSelection()
sce.subset<-devianceFeatureSelection(sce.subset, assay="counts",
                                     fam="poisson",sorted=TRUE)

pdf('plots/rank_vs_poissonDeviance_sceSubset.pdf')
plot(rowData(sce.subset)$poisson_deviance, type="l",
     xlab="ranked genes",
     ylab="poisson deviance",
     main="Feature Selection with Deviance",
     ylim=c(0,500000))
abline(v=3000,col='red',lty='dashed')
dev.off()
hdg<-rownames(counts(sce.subset))[1:3000]

#Run nullResiduals to calculate Pearson residuals from poisson model
sce.subset<-nullResiduals(sce.subset, assay = "counts",
                          fam = "poisson", type = "pearson")

#PCA
set.seed(2014)
sce.subset <- runPCA(sce.subset,exprs_values="poisson_pearson_residuals",
                     subset_row=hdg,ncomponents=50)
reducedDimNames(sce.subset)

set.seed(202014)
sce.subset <- runUMAP(sce.subset, dimred="PCA")

#Verify length
ncol(reducedDim(sce.subset, "PCA"))
sce.subset

pdf('plots/umap_mouse_sceSubset.pdf')
plotUMAP(sce.subset,colour_by='cellType',text_by='cellType')
dev.off()

# load NMF pattern
load("nmf_final.rda")

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Translate from one species to the other using the orthology
names <- orthology[orthology$Column3 %in% rownames(sce.subset),]

names <- names[match(rownames(sce.subset), names$Column3),]

setdiff(names$Column3, rownames(sce.subset))

rownames(sce.subset) <- names$Column1

spe <- sce.subset

# Project Erik's human snRNAseq NMF patterns onto mouse snRNAseq

set.seed(1029)
i<-intersect(rownames(spe),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]
proj<-project(logcounts(spe2),loadings,L1=0)

###scale projected weights
h <- (t(proj) / x@d)

###add to colData
colData(spe)<-cbind(colData(spe),h)

data<-as.data.frame(spe$cellType)
colnames(data)<-'cellType'
onehot_cellType <-  dcast(data = data, rownames(data) ~ cellType, length)
rownames(onehot_cellType)<-onehot_cellType[,1]
onehot_cellType[,1]<-as.numeric(onehot_cellType[,1])
onehot_cellType<-onehot_cellType[order(onehot_cellType[,1],decreasing=F),]
onehot_cellType[,1]<-NULL

correlation_celltype_projection <- as.matrix(cor(onehot_cellType, h))

colnames(correlation_celltype_projection) <- colnames(x@w)

# Find NAs
which(is.na(correlation_celltype_projection), arr.ind=TRUE)

#       row col
#GC.1     1   2
#GC.2     2   2
#CA4      3   2
#CA3.1    4   2
#CA3.2    5   2
#CA2      6   2
#CA1      7   2
#PS.1     8   2
#PS.2     9   2
#Sub     10   2
#L2/3    11   2
#L5/Po   12   2
#L6/6b   13   2
#GABA.1  14   2
#GABA.2  15   2
#GABA.3  16   2
#GABA.4  17   2
#GABA.5  18   2

# All from column 2 (NMF2) so remove that column

correlation_celltype_projection2 <- correlation_celltype_projection[,-2]

pdf(file = here::here("plots","NMF_human_to_mouse_heatmap.pdf"),
    width = 8.5, height = 9)

Heatmap(t(correlation_celltype_projection2),
    name = "corr",
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6))

dev.off()

sce.subset <- spe

saveRDS(sce.subset,
    file = here::here("sce_mouse_ECS.rds"))

# Plot UMAP and boxplot for nmf patters

###some duplicate barcodes, just paste sample to uniquify
colnames(sce.subset)<-paste(sce.subset$Sample,colnames(sce.subset),sep='_')
###group GCs
sce.subset$cellType<-factor(ifelse(sce.subset$cellType %in% c('GC.1','GC.2'),'GC',as.character(sce.subset$cellType)))

###plot boxplots for nmf patterns

pdf(file = here::here("plots","NMF_human_to_mouse_boxplots.pdf"))

ggcells(sce.subset,
        mapping=aes(x=cellType, y=V14,fill=condition)) +
        geom_boxplot(outlier.size = 0.05) +
    ylab("nmf14") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

ggcells(sce.subset,
        mapping=aes(x=cellType, y=V20,fill=condition)) +
        geom_boxplot(outlier.size = 0.05) +
    ylab("nmf20") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

ggcells(sce.subset,
        mapping=aes(x=cellType, y=V91,fill=condition)) +
        geom_boxplot(outlier.size = 0.05) +
    ylab("nmf91") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

ggcells(sce.subset,
        mapping=aes(x=cellType, y=V26,fill=condition)) +
        geom_boxplot(outlier.size = 0.05) +
    ylab("nmf26") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

dev.off()

# Plot UMAP

pdf(file = here::here("plots","NMF_human_to_mouse_UMAP.pdf"))

plotUMAP(sce.subset,colour_by='V14',text_by='cellType',point_size=0.1)+
    labs(color = "nmf14\nweight") +
scale_color_gradient(low = "grey", high = "black")

plotUMAP(sce.subset,colour_by='V20',text_by='cellType',point_size=0.1)+
    labs(color = "nmf20\nweight") +
scale_color_gradient(low = "grey", high = "black")

plotUMAP(sce.subset,colour_by='V91',text_by='cellType',point_size=0.1)+
    labs(color = "nmf91\nweight") +
scale_color_gradient(low = "grey", high = "black")

dev.off()

##############################################################################################

# Load macaque Hao scRNAseq dataset

sce <- readRDS(here::here("QCed_sce_hao.rds"))

# Run basic PCA and UMAP dimensionality reduction for plotting purposes

# Feature selection
dec <- modelGeneVar(sce, block = sce$sample_ID)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
# [1] 966
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 2861
head(top_hvgs)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 2399
head(top_hvgs)

# going with top.hvgs.fdr1

# Dimensionality reduction
set.seed(12345)
sce <- runPCA(sce, subset_row = top.hvgs.fdr1, ncomponents = 50)
set.seed(12345)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- c("UMAP1", "UMAP2")

# plot UMAP
pdf('plots/umap_macaque_sce.pdf')
plotUMAP(sce,colour_by='cellType',text_by='cellType')
dev.off()

# Project Erik's human snRNAseq NMF patterns onto macaque scRNAseq

set.seed(1029)
i<-intersect(rownames(sce),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
sce2<-sce[rownames(sce) %in% i,]
loadings<-loadings[match(rownames(sce2),rownames(loadings)),]
proj<-project(logcounts(sce2),loadings,L1=0)

###scale projected weights
h <- (t(proj) / x@d)

###add to colData
colData(sce)<-cbind(colData(sce),h)

data<-as.data.frame(sce$cellType)
colnames(data)<-'cellType'
onehot_cellType <-  dcast(data = data, rownames(data) ~ cellType, length)
rownames(onehot_cellType)<-onehot_cellType[,1]
onehot_cellType[,1]<-as.numeric(onehot_cellType[,1])
onehot_cellType<-onehot_cellType[order(onehot_cellType[,1],decreasing=F),]
onehot_cellType[,1]<-NULL

correlation_celltype_projection <- as.matrix(cor(onehot_cellType, h))

colnames(correlation_celltype_projection) <- colnames(x@w)

# Find NAs
which(is.na(correlation_celltype_projection), arr.ind=TRUE)

#            row col
#Astro_1       1  41
#Astro_2       2  41
#Astro_3       3  41
#Astro_4       4  41
#Astro_im1     5  41
#Astro_im2     6  41
#CR            7  41
#Endothelial   8  41
#Ependymal     9  41
#GABA         10  41
#GC_1         11  41
#GC_2         12  41
#GC_3         13  41
#GC_im        14  41
#IPC_1        15  41
#IPC_2        16  41
#Microglia_1  17  41
#Microglia_2  18  41
#MOL          19  41
#NB           20  41
#NFOL         21  41
#OPC_1        22  41
#OPC_2        23  41
#PhgMG_cRG    24  41
#PhgMG_Pyr    25  41
#PhgMG_RGL    26  41
#Pre_OPC      27  41
#PVM          28  41
#Pyr_1        29  41
#Pyr_2        30  41
#RGL_1        31  41
#RGL_2        32  41
#Unk          33  41
#VLMC         34  41

# All from column 41 (NMF41) so remove that column

correlation_celltype_projection2 <- correlation_celltype_projection[,-41]

correlation_celltype_projection2 <- t(correlation_celltype_projection2)

# Order columns to sort of match Eriks mouse celltypes ordered columns

list <- c("GC_1", "GC_2", "GC_3", "GC_im", "Pyr_1", "Pyr_2", "CR", "GABA",
    "NB", "IPC_1", "IPC_2", "RGL_1", "RGL_2", "Ependymal", "Astro_1",
    "Astro_2", "Astro_3", "Astro_4", "Astro_im1", "Astro_im2",
    "MOL", "NFOL", "OPC_1", "OPC_2", "Pre_OPC", "Microglia_1", "Microglia_2",
    "PhgMG_Pyr", "PhgMG_cRG", "PhgMG_RGL", "PVM", "Endothelial", "VLMC", "Unk")

correlation_celltype_projection3 <- correlation_celltype_projection2[,list]

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file = here::here("plots","NMF_human_to_macaque_heatmap.pdf"),
    width = 8.5, height = 9)

Heatmap(correlation_celltype_projection3,
    name = "corr",
    cluster_columns = FALSE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 6))

dev.off()

saveRDS(sce,
    file = here::here("sce_hao_macaque_humanNMFprojections.rds"))

###plot boxplots for nmf patterns

pdf(file = here::here("plots","NMF_human_to_macaque_boxplots.pdf"))

ggcells(sce,
        mapping=aes(x=cellType, y=V26,fill=cellType)) +
        geom_boxplot(outlier.size = 0.05)+
    ylab("nmf26 weights") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

dev.off()

###plot boxplots for nmf patterns 26, 5, and 14

nmf_patterns <- as.data.frame(colData(sce)[, c("cellType", "V26", "V5", "V14")])

long_data <- gather(nmf_patterns, key = "Variable", value = "Value", -cellType)

Variable_order <- c("V26", "V5", "V14")
long_data$Variable <- factor(long_data$Variable, levels = Variable_order)

colors <- c("V26" = "magenta", "V5" = "yellow", "V14" = "green")

# If one needs to remove some cell types and focus on a specific collection

long_data <- long_data[long_data$cellType == "GC_1" |
        long_data$cellType == "GC_im" |
        long_data$cellType == "NB" |
        long_data$cellType == "GC_2" |
        long_data$cellType == "RGL_1" |
        long_data$cellType == "RGL_2" |
        long_data$cellType == "IPC_1" |
        long_data$cellType == "IPC_2" |
        long_data$cellType == "GC_3"
        ,]

# Factor cellType to get propoer order on x-axis
long_data$cellType <- factor(long_data$cellType, levels = c("RGL_1", "RGL_2",
"IPC_1", "IPC_2", "NB", "GC_im", "GC_3", "GC_1", "GC_2"))

pdf(file = here::here("plots", "NMF","boxplots_select_nmf26514_macaque_new.pdf"))

# Create a boxplot using ggplot
ggplot(long_data, aes(x = cellType, y = Value, fill = Variable)) +
    geom_boxplot(color = "black", position = position_dodge(), outlier.size = 0.005) +
    labs(y = "weights") +
    coord_cartesian(ylim =  c(0, 2e-04)) +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

dev.off()

# Plot UMAP and rasterize for easy changes in Illustrator

pdf(file = here::here("plots","NMF_human_to_macaque_UMAP.pdf"))

ggrastr::rasterize(
plotUMAP(sce,colour_by='V26',text_by='cellType',point_size=0.1)+
    labs(color = "nmf26\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "magenta") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

ggrastr::rasterize(
plotUMAP(sce,colour_by='V5',text_by='cellType',point_size=0.1)+
    labs(color = "nmf5\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "yellow") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

ggrastr::rasterize(
plotUMAP(sce,colour_by='V14',text_by='cellType',point_size=0.1)+
    labs(color = "nmf14\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "green") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

dev.off()

hex_26 <- make_hexbin(sce, nbins = 100,
                   dimension_reduction = "UMAP", use_dims=c(1,2))

pdf(file = here::here("plots","NMF_human_to_macaque_schex_UMAP.pdf"))

plot_hexbin_meta(hex_26, col="V26", action="median")+
    labs(fill = "nmf26\nweight", x='UMAP1',y='UMAP2',title='nmf26 (immature granule cells)') +
    scale_fill_gradient(low = "grey", high = "black") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

plot_hexbin_meta(hex_26, col="V26", action="median")+
    labs(fill = "nmf26\nweight", x='UMAP1',y='UMAP2',title='nmf26 (immature granule cells)') +
    scale_fill_gradient(low = "grey", high = "magenta") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

plot_hexbin_meta(hex_26, col="V5", action="median")+
    labs(fill = "nm5\nweight", x='UMAP1',y='UMAP2',title='nmf5') +
    scale_fill_gradient(low = "grey", high = "yellow") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

plot_hexbin_meta(hex_26, col="V14", action="median")+
    labs(fill = "nm14\nweight", x='UMAP1',y='UMAP2',title='nmf14') +
    scale_fill_gradient(low = "grey", high = "green") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

dev.off()

# Try to get the color code consistent with UMAP
plotUMAPRGB<-function(sce, vars, ...) {
    plt_df <- data.frame(colData(sce), reducedDims(sce)$UMAP)

    if (any(!vars %in% names(plt_df))) {
        stop("One or more variables not found in the data.")
    }

    if (length(vars) > 4) {
        stop("A maximum of 4 variables is allowed.")
    }

    for (var in vars) {
        plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
    }

    num_vars <- length(vars)

    # Initialize channels based on the number of variables:
    # Magenta (R and B channels)
    # Yellow (R and G channels)
    # Green (G channel)
    # Blue (B channel)

    if (num_vars >= 1) {
        plt_df$R <- plt_df[[vars[1]]] # Part of Magenta and Yellow
        plt_df$B <- plt_df[[vars[1]]] # Part of Magenta
    }

    if (num_vars >= 2) {
        plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[vars[2]]] # Yellow component
        plt_df$G <- plt_df[[vars[2]]] # Green
    }

    if (num_vars >= 3) {
        plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[vars[3]]] # Green component
    }

    if (num_vars == 4) {
        plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[vars[4]]] # Blue component

    }

    sce$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
    plotUMAP(sce, colour_by = "RGB",...)+scale_colour_identity()
}

pdf(file = here::here("plots", "UMAP_nmf_26514gradient_macaque.pdf"))

ggrastr::rasterize(
plotUMAPRGB(sce,vars=c("V26", "V5", "V14"), text_by='cellType', point_size=0.1)
)

dev.off()

######################################################################################

# Hochgerner mouse postnatal neurogenesis scRNA-seq data

sce_hoch <- readRDS(here::here("QCed_sce_hochgerner.rds"))

# Run basic PCA and UMAP dimensionality reduction for plotting purposes

# Feature selection
dec <- modelGeneVar(sce_hoch, block = sce_hoch$strain)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
# [1] 1023
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 1275
head(top_hvgs)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 874
head(top_hvgs)

# going with top_hvgs

# Dimensionality reduction
set.seed(12345)
sce_hoch <- runPCA(sce_hoch, subset_row = top_hvgs, ncomponents = 50)
set.seed(12345)
sce_hoch <- runUMAP(sce_hoch, dimred = "PCA")
colnames(reducedDim(sce_hoch, "UMAP")) <- c("UMAP1", "UMAP2")

# plot UMAP
pdf('plots/umap_hochgerner_sce.pdf')
plotUMAP(sce_hoch,colour_by='cellType',text_by='cellType')
dev.off()

# Use table from Jax labs for mouse human orthology
#https://www.informatics.jax.org/downloads/reports/index.html#homology
orthology <- read.csv(file = here("processed-data","gene_set_enrichment",
    "human_mouse_orthologs.csv"))

# Translate from one species to the other using the orthology
names <- orthology[orthology$Column3 %in% rownames(sce_hoch),]

names <- names[match(rownames(sce_hoch), names$Column3),]

setdiff(names$Column3, rownames(sce_hoch))

#Replace NAs with string
names[] <- lapply(names, function(x) ifelse(is.na(x), "NA", x))

rownames(sce_hoch) <- names$Column1

# Project Erik's human snRNAseq NMF patterns onto mouse snRNAseq

set.seed(1029)
i<-intersect(rownames(sce_hoch),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
sce_hoch2<-sce_hoch[rownames(sce_hoch) %in% i,]
loadings<-loadings[match(rownames(sce_hoch2),rownames(loadings)),]
proj<-project(logcounts(sce_hoch2),loadings,L1=0)

###scale projected weights
h <- (t(proj) / x@d)

###add to colData
colData(sce_hoch)<-cbind(colData(sce_hoch),h)

data<-as.data.frame(sce_hoch$cellType)
colnames(data)<-'cellType'
onehot_cellType <-  dcast(data = data, rownames(data) ~ cellType, length)
rownames(onehot_cellType)<-onehot_cellType[,1]
onehot_cellType[,1]<-as.numeric(onehot_cellType[,1])
onehot_cellType<-onehot_cellType[order(onehot_cellType[,1],decreasing=F),]
onehot_cellType[,1]<-NULL

correlation_celltype_projection <- as.matrix(cor(onehot_cellType, h))

colnames(correlation_celltype_projection) <- colnames(x@w)

# Find NAs
which(is.na(correlation_celltype_projection), arr.ind=TRUE)

#                     row col
#Astro-adult            1   2
#Astro-juv              2   2
#CA3-Pyr                3   2
#Cajal-Retzius          4   2
#Endothelial            5   2
#Ependymal              6   2
#GABA                   7   2
#GC-adult               8   2
#GC-juv                 9   2
#Immature-Astro        10   2
#Immature-GABA         11   2
#Immature-GC           12   2
#Immature-Pyr          13   2
#MiCajal-Retziusoglia  14   2
#MOL                   15   2
#Neuroblast            16   2
#NFOL                  17   2
#nIPC                  18   2
#nIPC-perin            19   2
#OPC                   20   2
#PVM                   21   2
#RGL                   22   2
#RGL_young             23   2
#VLMC                  24   2
#Astro-adult            1  32
#Astro-juv              2  32
#CA3-Pyr                3  32
#Cajal-Retzius          4  32
#Endothelial            5  32
#Ependymal              6  32
#GABA                   7  32
#GC-adult               8  32
#GC-juv                 9  32
#Immature-Astro        10  32
#Immature-GABA         11  32
#Immature-GC           12  32
#Immature-Pyr          13  32
#MiCajal-Retziusoglia  14  32
#MOL                   15  32
#Neuroblast            16  32
#NFOL                  17  32
#nIPC                  18  32
#nIPC-perin            19  32
#OPC                   20  32
#PVM                   21  32
#RGL                   22  32
#RGL_young             23  32
#VLMC                  24  32
#Astro-adult            1  41
#Astro-juv              2  41
#CA3-Pyr                3  41
#Cajal-Retzius          4  41
#Endothelial            5  41
#Ependymal              6  41
#GABA                   7  41
#GC-adult               8  41
#GC-juv                 9  41
#Immature-Astro        10  41
#Immature-GABA         11  41
#Immature-GC           12  41
#Immature-Pyr          13  41
#MiCajal-Retziusoglia  14  41
#MOL                   15  41
#Neuroblast            16  41
#NFOL                  17  41
#nIPC                  18  41
#nIPC-perin            19  41
#OPC                   20  41
#PVM                   21  41
#RGL                   22  41
#RGL_young             23  41
#VLMC                  24  41

# All from column 41 (NMF41) so remove that column

correlation_celltype_projection2 <- correlation_celltype_projection[,-c(2,32,41)]

correlation_celltype_projection2 <- t(correlation_celltype_projection2)

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file = here::here("plots","NMF_human_to_hochgerner_heatmap.pdf"),
    width = 8.5, height = 9)

Heatmap(correlation_celltype_projection2,
    name = "corr",
    col_fun,
    row_names_gp = gpar(fontsize = 6))

dev.off()

saveRDS(sce_hoch,
    file = here::here("sce_hochgerner_mouse_humanNMFprojections.rds"))

###plot boxplots for nmf patterns

pdf(file = here::here("plots","NMF26_human_to_hochgerner_boxplots.pdf"))

ggcells(sce_hoch,
        mapping=aes(x=cellType, y=V26,fill=cellType)) +
        geom_boxplot(outlier.size = 0.05)+
    ylab("nmf26 weights") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

dev.off()

###plot boxplots for nmf patterns 26, 5, and 14

nmf_patterns <- as.data.frame(colData(sce_hoch)[, c("cellType", "V26", "V5", "V14")])

# Assuming your dataframe is named 'your_data'
# Reshape the data into long format
long_data <- gather(nmf_patterns, key = "Variable", value = "Value", -cellType)

Variable_order <- c("V26", "V5", "V14")
long_data$Variable <- factor(long_data$Variable, levels = Variable_order)

colors <- c("V26" = "magenta", "V5" = "yellow", "V14" = "green")

long_data <- long_data[long_data$cellType == "RGL" |
        long_data$cellType == "RGL_young" |
        long_data$cellType == "nIPC-perin" |
        long_data$cellType == "nIPC" |
        long_data$cellType == "Neuroblast" |
        long_data$cellType == "Immature-GC" |
        long_data$cellType == "GC-juv" |
        long_data$cellType == "GC-adult"
        ,]

# Factor cellType to get propoer order on x-axis
long_data$cellType <- factor(long_data$cellType, levels = c("RGL_young", "RGL",
"nIPC-perin", "nIPC", "Neuroblast", "Immature-GC", "GC-juv", "GC-adult"))

pdf(file = here::here("plots", "boxplots_select_nmf26514_hochgerner_new.pdf"))

# Create a boxplot using ggplot
ggplot(long_data, aes(x = cellType, y = Value, fill = Variable)) +
    geom_boxplot(color = "black", position = position_dodge(), outlier.size = 0.005) +
    labs(y = "weights") +
    coord_cartesian(ylim =  c(0, 1.5e-04)) +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

dev.off()

# Plot UMAP and rasterize for easy changes in Illustrator

pdf(file = here::here("plots","NMF_human_to_hochgerner_UMAP.pdf"))

ggrastr::rasterize(
plotUMAP(sce_hoch,colour_by='V26',text_by='cellType',point_size=0.1)+
    labs(color = "nmf26\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "magenta") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

ggrastr::rasterize(
plotUMAP(sce_hoch,colour_by='V5',text_by='cellType',point_size=0.1)+
    labs(color = "nmf5\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "yellow") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

ggrastr::rasterize(
plotUMAP(sce_hoch,colour_by='V14',text_by='cellType',point_size=0.1)+
    labs(color = "nmf14\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "green") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

dev.off()

hex <- make_hexbin(sce_hoch, nbins = 100,
                   dimension_reduction = "UMAP", use_dims=c(1,2))

pdf(file = here::here("plots","NMF_human_to_hochgerner_schex_UMAP.pdf"))

plot_hexbin_meta(hex, col="V26", action="median")+
    labs(fill = "nmf26\nweight", x='UMAP1',y='UMAP2',title='nmf26') +
    scale_fill_gradient(low = "grey", high = "black") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

plot_hexbin_meta(hex, col="V26", action="median")+
    labs(fill = "nmf26\nweight", x='UMAP1',y='UMAP2',title='nmf26') +
    scale_fill_gradient(low = "grey", high = "magenta") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

plot_hexbin_meta(hex, col="V5", action="median")+
    labs(fill = "nm5\nweight", x='UMAP1',y='UMAP2',title='nmf5') +
    scale_fill_gradient(low = "grey", high = "yellow") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

plot_hexbin_meta(hex, col="V14", action="median")+
    labs(fill = "nm14\nweight", x='UMAP1',y='UMAP2',title='nmf14') +
    scale_fill_gradient(low = "grey", high = "green") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

dev.off()

# Try to get the color code consistent with UMAP
plotUMAPRGB<-function(sce, vars, ...) {
    plt_df <- data.frame(colData(sce), reducedDims(sce)$UMAP)

    if (any(!vars %in% names(plt_df))) {
        stop("One or more variables not found in the data.")
    }

    if (length(vars) > 4) {
        stop("A maximum of 4 variables is allowed.")
    }

    for (var in vars) {
        plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
    }

    num_vars <- length(vars)

    # Initialize channels based on the number of variables:
    # Magenta (R and B channels)
    # Yellow (R and G channels)
    # Green (G channel)
    # Blue (B channel)

    if (num_vars >= 1) {
        plt_df$R <- plt_df[[vars[1]]] # Part of Magenta and Yellow
        plt_df$B <- plt_df[[vars[1]]] # Part of Magenta
    }

    if (num_vars >= 2) {
        plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[vars[2]]] # Yellow component
        plt_df$G <- plt_df[[vars[2]]] # Green
    }

    if (num_vars >= 3) {
        plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[vars[3]]] # Green component
    }

    if (num_vars == 4) {
        plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[vars[4]]] # Blue component

    }

    sce$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
    plotUMAP(sce, colour_by = "RGB",...)+scale_colour_identity()
}

pdf(file = here::here("plots", "UMAP_nmf_26514gradient_hochgerner.pdf"))

ggrastr::rasterize(
plotUMAPRGB(sce_hoch,vars=c("V26", "V5", "V14"), text_by='cellType', point_size=0.1)
)

dev.off()
