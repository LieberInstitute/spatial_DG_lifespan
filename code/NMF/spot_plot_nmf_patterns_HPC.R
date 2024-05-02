#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
#########plotting functions
library(SpatialExperiment)
library(scater)
library(here)

source(here::here('code', 'NMF','plotVisium_rewrite.R'))
source(here::here('code', 'NMF','visiumRGB_functions.R'))
load(here::here('processed-data', 'NMF','spe_nelson_nmf_final.rda'))
##set up speb for plotting
brains <- unique(spe$brnum)

speb <- spe[, (colData(spe)$brnum == brains[[5]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

###Set palette
spatial.palette3<-c("#BEDDBA", "#eae8e4", "#E8BBC6","#A1BAD8")
names(spatial.palette3)<-levels(speb$broad.domain)

# if plotting a gene instead, add logcounts to colData(spe)
spe$PPFIA2 <- assays(spe)$logcount["PPFIA2",]
spe$PRKCG <- assays(spe)$logcount["PRKCG",]
spe$APOC1 <- assays(spe)$logcount["APOC1",]
spe$SFRP2 <- assays(spe)$logcount["SFRP2",]

pdf(here::here('plots','spot_plots_Br2743_SFRP2_.pdf'))

###make single-variable plot with broad domain ground
ggrastr::rasterize(
plotVisium(
    speb,
    spots = TRUE,
    fill = 'SFRP2',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = spatial.palette,
    image=F,
    values=spatial.palette3)+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
)

dev.off()

pdf(here::here('plots', 'NMF','spot_plots_Br6471_nmf26_5_14.pdf'))

###make multi-variable plot (var 1=magenta, var 2=yellow, var 3=green, var 4=blue)
plotVisiumRGB(speb, pink = 'nmf26', yellow = 'nmf5', green = 'nmf14',
image=F,highlight='neuron_cell_body',values=c('gray50','black'))

dev.off()
