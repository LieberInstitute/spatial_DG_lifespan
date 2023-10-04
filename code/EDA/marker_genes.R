###############################
# spatial_DG_lifespan project
# Plotting marker genes
# Anthony Ramnauth, Apr 21 2022
###############################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(ggnewscale)
    library(spatialLIBD)
    library(ggspavis)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "marker_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

# order spe observations according to age
spe <- spe[, order(spe$age)]

# Find marker genes
human_markers <-
    c(
        "SNAP25",
        "SLC17A7",
        "GAD1",
        "GAD2",
        "MBP",
        "MOBP",
        "RELN",
        "AQP4",
        "CCK",
        "HPCAL1",
        "PROX1",
        "NECAB1",
        "MPPED1",
        "SLC17A6",
        "TNNT2",
        "CALB1",
        "GABRQ",
        "THOC3",
        "KCNJ4",
        "SFRP2",
        "APLNR",
        "CD44",
        "TMEM155",
        "SYT13",
        "SLC25A22",
        "SLIT1",
        "SYT4",
        "APOC1",
        "WIF1",
        "MTRNR2L10",
        "TCF7L2",
        "RORA",
        "LEF1",
        "ZMAT4"
    )

# Locate the marker genes
human_markers_search <- rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_name)]

# Plot marker genes on tissue
for (i in human_markers_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "logcounts",
        minCount = 0,
        viridis = FALSE,
        point_size = 1.5
    )
}

# Find neurogenesis marker genes
neurogenesis_markers <-
    c(
        "DCX",
        "SOX11",
        "NEUROD1",
        "NEUROD2",
        "SOX2",
        "NEUROG2",
        "MCM2",
        "MKI67",
        "ELAVL2",
        "TACC2",
        "PALLD",
        "NNAT",
        "DPYSL5",
        "DPYSL3",
        "KLF7",
        "BHLHE22",
        "POSTN",
        "BTBD3",
        "CD24",
        "NPDC1",
        "DNM1",
        "NCDN",
        "CAMK2B",
        "NPTX1",
        "RFX3",
        "MALAT1",
        "GFAP",
        "APC"
    )

# Locate the neurogenesis genes
neurogenesis_markers_search <- rowData(spe)$gene_search[match(neurogenesis_markers, rowData(spe)$gene_name)]

# Plot neurogenesis marker genes on tissue
for (i in neurogenesis_markers_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "logcounts",
        minCount = 0,
        viridis = FALSE,
        point_size = 1.5
    )
}

## Trying out custom function for "fluorescent-like" multiple expression labeling of histology

# For each gene of interest you need to create a column in colData(spe) of their logcounts
# You can use Ensembl IDs or convert the rownames(spe) to

rownames(spe) <- rowData(spe)$gene_name

colData(spe)$GAD2 <- assays(spe)$logcounts["GAD2",]
colData(spe)$KIT <- assays(spe)$logcounts["KIT",]
colData(spe)$CALB1 <- assays(spe)$logcounts["CALB1",]
colData(spe)$TRHDE <- assays(spe)$logcounts["TRHDE",]
colData(spe)$AATK <- assays(spe)$logcounts["AATK",]
colData(spe)$SERPINI1 <- assays(spe)$logcounts["SERPINI1",]

# Add column in colData() for neuron
spe$neuron <- FALSE
spe$neuron[spe$bayesSpace_harmony_10 == "4" |
        spe$bayesSpace_harmony_10 == "7"|
        spe$bayesSpace_harmony_10 == "9"] <- TRUE

plotVisiumRGB <- function(spe, vars, ...) {
    plt_df <- data.frame(colData(spe), spatialCoords(spe))

    if (any(!vars %in% names(plt_df))) {
        stop("One or more variables not found in the data.")
    }

    for (var in vars) {
        plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
    }

    num_vars <- length(vars)

    # Initialize RGB channels based on number of variables
    if (num_vars >= 3) {
        plt_df$R <- plt_df[[vars[1]]]
        plt_df$G <- plt_df[[vars[2]]]
        plt_df$B <- plt_df[[vars[3]]]
    } else {
        plt_df$R <- plt_df[[vars[1]]]
        plt_df$B <- plt_df[[vars[2]]]
        plt_df$G <- rep(0, nrow(plt_df))
    }

    if (num_vars >= 4) {
        W <- plt_df[[vars[4]]]
        plt_df$R <- plt_df$R + W * (1 - plt_df$R)
        plt_df$G <- plt_df$G + W * (1 - plt_df$G)
        plt_df$B <- plt_df$B + W * (1 - plt_df$B)
    }

    if (num_vars >= 5) {
        C <- plt_df[[vars[5]]]
        plt_df$G <- plt_df$G + C * (1 - plt_df$G)
        plt_df$B <- plt_df$B + C * (1 - plt_df$B)
    }

    if (num_vars == 6) {
        M <- plt_df[[vars[6]]]
        plt_df$R <- plt_df$R + M * (1 - plt_df$R)
        plt_df$B <- plt_df$B + M * (1 - plt_df$B)
    }

    spe$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
    plotVisium(spe, fill = "RGB", ...)+scale_fill_identity()
}

# Plot new spot plots with marker genes

plotVisiumRGB(spe, c("PPFIA2", "SPOCK1", "MET", "BACE1", "PTPRT", "GAD2"), highlight = "neuron",
    image = FALSE, facets = "age")

pdf(file = here::here("plots", "marker_genes", "TEST_coexpression.pdf"), width = 12, height = 12)

plotVisiumRGB(spe, c("GAD2", "CALB1", "KIT"), highlight = "neuron",
    image = FALSE, facets = "age") +
    theme(legend.position="none")

dev.off()

