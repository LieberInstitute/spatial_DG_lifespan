################################
# spatial_HPC project
# Pseudobulk SPE by capture area
# Anthony Ramnauth, Ded 04 2023
################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(ggplot2)
    library(ggnewscale)
    library(scater)
    library(scran)
    library(dplyr)
    library(tidyr)
    library(ggsignif)
    library(CoGAPS)
    library(RcppML)
    library(projectR)
    library(reshape2)
    library(circlize)
    library(ComplexHeatmap)
    library(schex)
})

# Load SCE
load(file=here::here('sce_nmf_final.rda'))

features1 <- c("GAD1", "LAMP5", "TAC1", "CCK", "PDYN")

pdf(file = here::here("violinplot_GAD1_LAMP5_TAC1_CCK_PDYN.pdf"), width = 12, height = 12)

plotExpression(sce, x = "fine.type", features = features1, colour_by = "fine.type", ncol = 1) +
    geom_boxplot(width = 0.1) +
    theme(text = element_text(size = 20), axis.text = element_text(size = 14),
        legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic"))

dev.off()

###################################################################################
# Explore nmfs within this dataset

data<-as.data.frame(sce$fine.type)
colnames(data)<-'cellType'
onehot_cellType <-  dcast(data = data, rownames(data) ~ cellType, length)
rownames(onehot_cellType)<-onehot_cellType[,1]
onehot_cellType[,1]<-as.numeric(onehot_cellType[,1])
onehot_cellType<-onehot_cellType[order(onehot_cellType[,1],decreasing=F),]
onehot_cellType[,1]<-NULL

correlation_celltype_projection <- as.matrix(cor(onehot_cellType, as.matrix(colData(sce)[,32:131])))
correlation_celltype_projection <- t(correlation_celltype_projection)

# Order columns to sort of match Eriks mouse celltypes ordered columns
# If I have the time to polish
#list <- c()

#correlation_celltype_projection3 <- correlation_celltype_projection2[,list]

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file = here::here("plots","NMF_HPC_chromium_snRNA_heatmap.pdf"),
    width = 9, height = 9)

Heatmap(correlation_celltype_projection,
    name = "corr",
    col = col_fun,
    row_names_gp = gpar(fontsize = 6))

dev.off()

# isolate the GCs
sce_gc <- sce[, sce$cell.type == "GC"]

# Plot UMAP for nmf26
pdf(file = here::here("plots", "UMAP_nmf26_HPC_Chromium.pdf"))

ggrastr::rasterize(
plotUMAP(sce,colour_by='nmf26',text_by='fine.type',point_size=0.1)+
    labs(color = "nmf26\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "grey", high = "black") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

hex_26 <- make_hexbin(sce, nbins = 100,
                   dimension_reduction = "UMAP", use_dims=c(1,2))

plot_hexbin_meta(hex_26, col="nmf26", action="median")+
    labs(fill = "nmf26\nweight", x='UMAP1',y='UMAP2',title='nmf26 (immature granule cells)') +
    scale_fill_gradient(low = "grey", high = "black") +
    theme(axis.ticks=element_blank(),axis.text=element_blank())

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce$cell.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "cell.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce$nmf26,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf26") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_gc, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce_gc$fine.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "fine.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_gc, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce_gc$nmf26,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf26") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

dev.off()

# Plot UMAP for nmf20
pdf(file = here::here("plots", "UMAP_nmf20_HPC_Chromium.pdf"))

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce$cell.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "cell.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce$nmf20,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf20") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_gc, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce_gc$fine.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "fine.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_gc, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce_gc$nmf20,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf20") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

dev.off()

# isolate the NeuroImmune cluster
sce_immune <- sce[, sce$cell.type == "Micro/Macro/T"]

# Plot UMAP for nmf90
pdf(file = here::here("plots", "UMAP_nmf90_HPC_Chromium.pdf"))

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce$cell.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "cell.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce$nmf90,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf90") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_immune, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce_immune$fine.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "fine.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_immune, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce_immune$nmf90,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf90") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

dev.off()

# isolate the GABA cluster
sce_GABA <- sce[, sce$cell.type == "GABA"]

# Plot UMAP for nmf73
pdf(file = here::here("plots", "UMAP_nmf73_HPC_Chromium.pdf"))

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce$cell.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "cell.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce$nmf73,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf73") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_GABA, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = factor(sce_GABA$fine.type),
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "fine.type") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce_GABA, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce_GABA$nmf73,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf73") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()

dev.off()

# Plot UMAP for gene markers from nmf patterns

plotReducedDim(
    sce,
    dimred = "UMAP",
    ncomponents = 2,
    colour_by = "POSTN",
    point_size = 0.1,
    point_alpha = 0.5
) +
    scale_color_gradient(low = "grey", high = "black") +
    labs(color = "POSTN") +
    theme_bw() +
    theme_classic()

plotReducedDim(
    sce,
    dimred = "UMAP",
    ncomponents = 2,
    point_size = 0.1,
    colour_by = "COL25A1",
    point_alpha = 0.5
) +
    scale_color_gradient(low = "grey", high = "black") +
    labs(color = "COL25A1") +
    theme_bw() +
    theme_classic()

ggplot(data.frame(reducedDim(sce, "UMAP")),
    aes(
        x = UMAP1,
        y = UMAP2,
        color = sce$nmf69,
        alpha = 0.5
    )) +
    geom_point() +
    labs(color = "nmf90") +
    scale_color_gradient(low = "grey", high = "black") +
    theme_bw() +
    theme_classic()


###plot boxplots for nmf patterns

nmf_patterns <- as.data.frame(colData(sce)[, c("cell.type", "nmf14", "nmf20", "nmf26", "nmf91")])

# Assuming your dataframe is named 'your_data'
# Reshape the data into long format
long_data <- gather(nmf_patterns, key = "Variable", value = "Value", -cell.type)

colors <- c("nmf14" = "blue", "nmf20" = "darkgreen", "nmf26" = "darkgrey", "nmf91" = "lightgoldenrod")

pdf(file = here::here("plots", "boxplots_select_nmf_HPC_Chromium.pdf"))

# Create a boxplot using ggplot
ggplot(long_data, aes(x = cell.type, y = Value, fill = Variable)) +
    geom_boxplot(color = "black", position = position_dodge(), outlier.size = 0.05, width = 1) +
    labs(y = "weights") +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

dev.off()

# Boxplot of only nmf26

ggcells(sce,
        mapping=aes(x=fine.type, y=nmf26)) +
        geom_boxplot(outlier.size = 0.05, fill = "grey") +
    ylab("nmf26") +
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

# Make RGB UMAP for patterns 26, 5, and 14
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

pdf(file = here::here("plots", "NMF","UMAP_nmf_26514gradient.pdf"))

ggrastr::rasterize(
plotUMAPRGB(sce,vars=c("nmf26", "nmf5", "nmf14"), text_by='fine.type', point_size=0.1)
)

dev.off()

# Make heatmap of weights for genes of select nmf patterns

# Load marker genes for NMF patterns
load("patternmarks.rda")
# load NMF pattern
load("nmf_nelson_final.rda")

marks$PatternMarkers$nmf14
#[1] "GALNT13"    "C8orf34"    "ZMAT4"      "PTCHD1-AS"  "ACVR1C"     "LINC02248"  "PLEKHA2"    "NPY1R"      "AC113383.1"
#"PTCHD1"
#[11] "AC010967.1" "AHR"        "LINC01931"  "AC005323.2" "VWA5A"      "MYO3B"      "TGFA"       "COL27A1"    "SLITRK3"
#"AC009878.2"
#[21] "NECAB3"     "AC004817.5" "P3H2-AS1"   "AC087627.1" "AC117944.1" "SDCBP2"     "AL138799.4" "CENPI"      "LRRC53"
#"AC004470.1"
#[31] "AL359764.3" "AL359764.2" "AQP3"       "MPPED2-AS1" "AL359851.1" "LINC01121"  "AC103796.1" "OPN4"       "LINC01847"
#"AL662844.4"
#[41] "AC022915.2" "AL121929.2" "AC093843.2" "AC073548.2" "MYO18B"     "AL356580.1" "MIR4453HG"  "AC022196.1" "PRSS53"
#"AC009812.1"
#[51] "AC010291.1" "AC068481.1" "PHEX-AS1"   "AC022424.2" "AC023510.1" "AC004817.1" "AC004828.2" "AC017002.5" "AC092832.2"
#"LINC01707"
#[61] "AC005323.1" "Z82196.2"   "AC009478.1" "LINC01498"  "AC004828.1" "TRMT2B-AS1" "AC108863.2" "AC004817.4" "AC018643.1"
#"HTRA4"
#[71] "AC093158.1" "AC022272.1" "NOX3"       "AC018464.1" "ADGRG4"     "AL080284.1" "PGA3"       "AC138965.1"

marks$PatternMarkers$nmf20
# [1] "SORCS3"     "PTPRO"      "SPRY4-AS1"  "BDNF"       "MICAL2"     "RASGRF1"    "OTOGL"      "FZD3"       "ETV5"
#"RPS6KA3"
#[11] "HRH1"       "VGF"        "LINC02073"  "AC087280.2" "DLGAP4"     "RNF24"      "RPS6KA6"    "NPTX2"      "SLC6A17"
#"AC002377.1"
#[21] "PTPN5"      "GALNT9"     "PCSK1"      "FBXO16"     "INSYN2A"    "CYTOR"      "SPNS2"      "HSD17B12"   "SYT12"
#"ERCC1"
#[31] "CAMK1G"     "SPRY2"      "DUSP6"      "GCNA"       "ZDHHC5"     "RNF128"     "UBR7"       "FRMD8"      "PLPBP"
#"CBARP"
#[41] "C2orf92"    "DUSP14"     "AC112722.1" "SPATA2L"    "AC090125.1" "INHBA"      "MEG9"       "SPRED3"     "CCKBR"
#"AC016727.1"
#[51] "WDR55"      "DNAJB5"     "TAF5"       "AC068205.2" "STAC3"

marks$PatternMarkers$nmf26
#[1] "ALDH1A2"    "AC079793.1" "FIGN"       "TARID"      "POSTN"      "WIPF3"      "NREP"       "LINC01885"  "AL590807.1"
#"RHBDL3"
#[11] "CALB1"      "LINC00871"  "LINC01151"  "LONRF3"     "MCOLN3"     "PTGFR"      "EPHA5-AS1"  "KLHDC8A"    "ZDHHC22"
#"SMIM17"
#[21] "AC009975.1" "AJ006995.1" "TMEM244"    "BHLHE22"    "TDRD9"      "ZNF107"     "NEB"        "AKAIN1"     "ZNF682"
#"AC061961.1"
#[31] "PDYN"       "SLC9A5"     "LINC02338"  "SLC27A2"    "CCDC77"     "AC103855.2" "LINC02343"  "ADGRG2"     "LINC00457"
#"SLC1A6"
#[41] "TRIM35"     "AC018541.1" "SLC35F2"    "PABPC1L"    "TRHR"       "HLX-AS1"    "LINC01695"  "SGCG"       "CCDC169"
#"AL591368.1"
#[51] "AC093898.1" "SYT3"       "CLSPN"      "AC104248.1" "LINC01338"  "P2RX2"      "AC113346.2" "AC103855.3" "IL27RA"
#"NANP"
#[61] "AP001496.4" "AC105429.1" "AQP9"       "AC112721.2" "AC009118.3" "NPY2R"      "AL161891.1" "MYCN"       "DUSP9"
#"AC006197.2"
#[71] "LINC01315"  "AC104024.2" "PMFBP1"     "LINC02099"  "AC104699.1" "LINC02070"  "ZCCHC18"    "AC104407.1" "AC009230.1"
#"AP000248.1"
#[81] "AC092651.3" "LINC01789"  "FST"        "AC096669.1" "LINC02484"  "AC092100.1" "AC082651.3" "AL132633.1" "STC1"
#"AC016687.3"
#[91] "LINC02715"  "LINC02775"  "AC073048.1" "ZBTB20-AS3" "LINC00314"  "AL024497.1" "AL024497.2" "SMILR"

marks$PatternMarkers$nmf91
#[1] "SAMD4A"     "ZBTB16"     "PDK4"       "SIK2"       "FGFR1"      "ELOVL5"     "NR4A3"      "ESYT2"      "FAT1"
#"TLE1"
# [11] "JUND"       "HMGCS1"     "JUN"        "CTNNA1"     "MXI1"       "FOS"        "ARHGEF10"   "SPAG9"      "TFRC"
#"ZNF331"
# [21] "ATP1B3"     "DUSP16"     "ABHD3"      "NFE2L2"     "SREBF2"     "CCND2"      "IFRD1"      "ZFAND5"     "ZFAND6"
#"ARIH1"
# [31] "RHEB"       "RHOQ"       "RNF19A"     "NPC1"       "CNN3"       "SLC38A2"    "FDFT1"      "STAT3"      "MRTFA"
#"DNAJB6"
# [41] "MIRLET7BHG" "NR4A1"      "NARS2"      "NAMPT"      "BNIP3L"     "BCL6"       "UBE2B"      "HSPD1"      "MSMO1"
#"LSS"
# [51] "NR4A2"      "AGO2"       "IDI1"       "HNRNPU"     "CYCS"       "RANBP2"     "LONRF1"     "LHX2"       "SQSTM1"
#"MORF4L2"
# [61] "SKI"        "VPS37B"     "GEM"        "INSIG1"     "GABARAPL1"  "RAB11FIP3"  "CYSTM1"     "TSPYL2"     "SAMD4B"
#"CRY1"
# [71] "FADS2"      "KLF9"       "ZNF10"      "ACO2"       "TPM4"       "NUFIP2"     "PPP4R3A"    "GATAD2A"    "TOB1"
#"HSPA9"
# [81] "TMEM184B"   "PRMT9"      "CLK1"       "ANP32A"     "EGR1"       "LDLR"       "NUMA1"      "TXNL4A"     "MED14"
#"MAT2A"
# [91] "HSPE1"      "CYP51A1"    "TTC19"      "FOXK2"      "CCNL1"      "AZIN1"      "PITPNA"     "TMED10"     "JUNB"
#"BTG2"
#[101] "HSPA5"      "RAB21"      "HBP1"       "TOB2"       "SLC20A1"    "SYAP1"      "MPC1"       "TMEM41B"    "NUP98"
#"CRY2"
#[111] "SLC25A33"   "NAPA"       "MXD1"       "PPM1D"      "ARF4"       "SQLE"       "ALKBH5"     "PER1"       "FADS1"
#"ZNF141"
#[121] "EIF4H"      "SRSF3"      "ARFGAP3"    "EIF4A1"     "NDEL1"      "YTHDF3"     "POLR2A"     "USP36"      "SLC66A2"
#"SC5D"
#[131] "BCAS2"      "IRF2BP2"    "PPP2R1B"    "MVD"        "ARMCX3"     "BHLHE40"    "BRD2"       "WDR45B"     "XRCC6"
#"MOB4"
#[141] "RNF138"     "TXLNG"      "FDPS"       "DNAJB1"     "AKIRIN1"    "PIP5K1C"    "YPEL5"      "PER2"       "MAPK1IP1L"
#"ADNP2"
#[151] "SLC33A1"    "ETF1"       "UAP1"       "TOR1AIP1"   "AC005261.1" "RBBP7"      "HAT1"       "JOSD1"      "TSC22D3"
#"ZNF891"
#[161] "BIRC2"      "FAM210A"    "C1orf198"   "FOSL2"      "LMNA"       "KAT7"       "SLC25A25"   "SS18L1"     "ZNF791"
#"NFKBIZ"
#[171] "TRIM26"     "USP2"       "PHF10"      "FASN"       "HES1"       "MAPKAPK2"   "ACAT2"      "KLHL21"     "TIPARP"
#"TLE3"
#[181] "RIPK2"      "HNRNPF"     "STARD4"     "ZNF410"     "SLC19A2"    "SCML1"      "DDX27"      "CDC42SE1"   "FAM53C"
#"ERRFI1"
#[191] "STX5"       "STX17-AS1"  "CCNT1"      "MIDN"       "ARRDC3"     "SRA1"       "DDIT4"      "OPA3"       "ERVK3-1"
#"GTPBP1"
#[201] "CDK11B"     "IER2"       "LMBR1L"     "COQ10B"     "IFFO2"      "DESI1"      "ING3"       "GNL3L"      "ZBTB21"
#"MNT"
#[211] "CASP9"      "GNL3"       "SIDT2"      "FN3K"       "ISG20L2"    "PMVK"       "ANKRD37"    "INPP5K"     "SMNDC1"
#"CTH"
#[221] "HSD17B7"    "CCN1"       "NECTIN2"    "ERF"        "YBX3"       "PPP1R15B"   "TMEM104"    "C4orf33"    "TBC1D25"
#"ZBTB5"
#[231] "NR1H2"      "EHD1"       "MAPRE1"     "KLHL15"     "AC026471.4" "STK19"      "BCORL1"     "OSER1"      "EMD"
#"ZNF778"
#[241] "APBB3"      "PLA2G6"     "ZNF787"     "GMEB2"      "MKNK2"      "BAG3"       "SRF"        "GZF1"       "MIR22HG"
#"RIT1"
#[251] "CNTROB"     "XBP1"       "ZNF548"     "UNG"        "ATF3"       "POLM"       "RRNAD1"     "SPART-AS1"  "NSDHL"
#"FOSB"
#[261] "SNHG15"     "GADD45B"    "PPT2"       "CSRNP1"     "PIGA"       "RGS16"      "AC011379.1" "MAFF"       "AC104695.4"

# Filter out low weight genes for each nmf pattern
nmf_14_df <- as.matrix(x@w[marks$PatternMarkers$nmf14[1:20], c("nmf14")])

summary(nmf_14_df)
#       V1
# Min.   :0.0001322
# 1st Qu.:0.0003216
# Median :0.0004635
# Mean   :0.0005539
# 3rd Qu.:0.0006910
# Max.   :0.0011920

nmf_14_features <- rownames(nmf_14_df)[which(nmf_14_df > 0.0005539, arr.ind = TRUE)[, 1]]

nmf_20_df <- as.matrix(x@w[marks$PatternMarkers$nmf20[1:20], c("nmf20")])

summary(nmf_20_df)
#       V1
# Min.   :0.0002320
# 1st Qu.:0.0005571
# Median :0.0006532
# Mean   :0.0008117
# 3rd Qu.:0.0010649
# Max.   :0.0020123

nmf_20_features <- rownames(nmf_20_df)[which(nmf_20_df > 0.0008117, arr.ind = TRUE)[, 1]]

nmf_26_df <- as.matrix(x@w[marks$PatternMarkers$nmf26, c("nmf26")])

summary(nmf_26_df)
#       V1
# Min.   :0.0002587
# 1st Qu.:0.0003686
# Median :0.0005365
# Mean   :0.0005994
# 3rd Qu.:0.0008337
# Max.   :0.0012451

nmf_26_features <- rownames(nmf_26_df)[which(nmf_26_df > 2.080e-04, arr.ind = TRUE)[, 1]]

nmf_91_df <- as.matrix(x@w[marks$PatternMarkers$nmf91[1:20], c("nmf91")])

summary(nmf_91_df)
#       V1
# Min.   :0.0007326
# 1st Qu.:0.0011622
# Median :0.0018242
# Mean   :0.0019370
# 3rd Qu.:0.0025049
# Max.   :0.0047864

nmf_91_features <- rownames(nmf_91_df)[which(nmf_91_df > 0.0019370, arr.ind = TRUE)[, 1]]

nmf_features <- unique(c(nmf_14_features, nmf_20_features, nmf_26_features, nmf_91_features))

nmf_df <- x@w[c(nmf_26_features), c("nmf26")]

col_fun1 = colorRamp2(c(0, 0.0015), c("white", "black"))

pdf(file = here::here("plots","NMF26_select_chromium_snRNA_heatmap.pdf"), height = 2)

Heatmap(t(nmf_df),
    name = "nmf26\nweight",
    col = col_fun1,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    column_names_gp = gpar(fontface = "italic"))

dev.off()

# Make heatmap for nmf26, nmf5, and nmf10 top 5 protein coding marker genes each

nmf_26_genes <- row.names(as.matrix(x@w[marks$PatternMarkers$nmf26, c("nmf26")]))
nmf_5_genes <- row.names(as.matrix(x@w[marks$PatternMarkers$nmf5, c("nmf5")]))
nmf_14_genes <- row.names(as.matrix(x@w[marks$PatternMarkers$nmf14, c("nmf14")]))

# Manually check which are protein coding (checking with www.genecards.org)

gene_list <- c("ALDH1A2", "FIGN", "POSTN", "WIPF3", "NREP",
    "SEMA5A", "ADRA1A", "TLL1", "DPF3", "ADAMTS16",
    "GALNT13", "C8orf34", "ZMAT4", "ACVR1C", "PLEKHA2")


nmf_df1 <- x@w[c(gene_list), c("nmf26", "nmf5", "nmf14")]

col_fun2 = colorRamp2(c(0, 0.002), c("white", "black"))

pdf(file = here::here("plots","NMF26514_top5proteincoding_chromium_snRNA_heatmap.pdf"))

Heatmap(nmf_df1,
    name = "nmf\nweight",
    col = col_fun2,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontface = "italic"))

dev.off()

# plot all of the marker genes to check how they look
nmf_df2 <- x@w[unique(c(nmf_26_genes, nmf_5_genes, nmf_14_genes)), c("nmf26", "nmf5", "nmf14")]

col_fun2 = colorRamp2(c(0, 0.002), c("white", "black"))

pdf(file = here::here("plots","NMF26514_all_chromium_snRNA_heatmap.pdf"))

Heatmap(nmf_df2,
    name = "nmf\nweight",
    col = col_fun2,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    row_names_gp = gpar(fontface = "italic"))

dev.off()

# Lets plot the top 10 makers per nmf pattern
gene_list2 <- c(head(marks$PatternMarkers$nmf26,10), head(marks$PatternMarkers$nmf5,10),
    head(marks$PatternMarkers$nmf14,10))

nmf_df2 <- x@w[c(gene_list2), c("nmf26", "nmf5", "nmf14")]

col_fun3 = colorRamp2(c(0, 0.002), c("white", "black"))

pdf(file = here::here("plots","NMF26514_top10_chromium_snRNA_heatmap_test.pdf"), width = 3)

Heatmap(nmf_df2,
    name = "nmf\nweight",
    col = col_fun3,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontface = "italic"))

dev.off()

# dotplot of canonical neurogenesis markers form GSEA

sce_GC <- sce[,sce$cell.type == "GC"]

# features taken from  GSEA
features <- c("ALDH1A2", "CNGB1", "PROX1", "MCOLN3", "PDLIM5", "LPAR3", "FAT4", "DOCK10", "NEUROD1", "SLIT3",
    "PLXNA4", "BCL11B", "NEUROD2", "PARD3", "TSPAN2", "USH2A", "FLRT3", "PRDM8", "TRPC5", "NPTX1",
    "NRP1", "GABRA5", "LAMB1", "PTPRO", "NTNG2", "TRPC4",
    "LRTM2",
    "TIAM1", "NLGN4X", "MYCN", "VAX2", "EPHA5", "MT1X", "GDF7", "ENC1",
    "DCC",
    "NREP", "APCDD1", "PTK2B", "RASAL1", "NYAP2", "CAMK2A", "KCNA1", "SDK1",
    "ATP8B1",
    "HHIP", "RGS14", "MYT1L", "ISLR2", "DTNB", "WASF1", "SDC2", "ARHGAP44",
    "LHX2",
    "LHFPL5", "TNIK", "NEGR1", "FEZF2", "ZHX2", "CAMK2B", "EHD1", "DISC1", "KALRN", "ADCY1",
    "MMD", "IL1RAPL1", "NRP2", "TNR", "ARK2C", "FGF5", "BTBD3", "NEDD4", "SHANK1",
    "SYT17",
    "ITPKA", "LRP12", "DCLK1", "PRDM1", "TENM1", "ITGA4", "HDAC1", "PLPPR4", "CDH4",
    "BOC",
    "NOTCH2", "NTRK3", "GABRB2", "NOG", "RARB", "RNF112", "PCP4", "OPRM1", "SEZ6", "GPRIN1",
    "LGI1",
    "PAK6", "HCN1", "SYT3", "GABRB1", "BICDL1", "CELSR2", "PPP3CA", "PARD6B",
    "PAQR3",
    "CDK5R2", "CDON", "SYNGAP1", "EFHC2", "CDK5R1", "FZD1")

# After reviewing, slim down to those with high z-scores for GC.3

features <- c("LPAR3", "DOCK10", "NEUROD2", "TSPAN2", "PRDM8", "TRPC4", "NTNG2", "GDF7",
    "APCDD1", "HHIP", "FGF5", "LRP12", "BOC", "NOG", "RARB", "OPRM1", "SEZ6", "FZD1")

features <- features[! features %in% setdiff(features, rownames(sce_GC))]

pdf(file = here::here("plots","Neurogenesis_markers_GCs.pdf"), height = 4)

plotDots(sce_GC, group = "fine.type", features = features, exprs_values = "logcounts", color = c("blue", "white", "red"),
    block = "brnum", center = TRUE, scale = TRUE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic")) +
    coord_flip()

dev.off()
