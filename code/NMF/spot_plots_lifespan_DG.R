#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
#########plotting functions
library(SpatialExperiment)
library(scater)
library(dplyr)
library(ggh4x)
library(scater)
library(here)

source(here::here('code', 'NMF', 'plotVisium_rewrite.R'))
source(here::here('code', 'NMF', 'visiumRGB_functions.R'))
spe <- readRDS(here::here('processed-data', 'NMF', 'spe_lifespan_DG_NMF_projections.rds'))

##set up broad domain for spot plots for plotting
spe$broad.domain <-factor(ifelse(spe$bayesSpace_harmony_10 %in% c(4,7,9),'Neuron',
                  ifelse(spe$bayesSpace_harmony_10 %in% c(1,2,5,6,8),'Neuropil',
                         ifelse(spe$bayesSpace_harmony_10 %in% c(10),'WM',
                          'Vasc_CSF'))))

###Set palette
spatial.palette3<-c("#BEDDBA", "#eae8e4", "#A1BAD8", "#E8BBC6")
names(spatial.palette3)<-levels(spe$broad.domain)

##set up speb for plotting individual capture areas

# Ensure that each barcode is unique use spe$key for colnames
colnames(spe) <- spe$key

# Add variable of age_bin to colData(spe)
age_df <- data.frame(spe$key, spe$sample_id, spe$age)
age_df <- age_df %>%
    mutate(age_bin = case_when(
        between(spe.age, 0, 3) ~ "Infant",
        between(spe.age, 13, 19) ~ "Teen",
        between(spe.age, 20, 50) ~ "Adult",
        between(spe.age, 60, 100) ~ "Elderly"
    ))

stopifnot(age_df$spe.key == spe$key)

colData(spe)$age_bin <- factor(age_df$age_bin, levels = c("Infant", "Teen", "Adult", "Elderly"))

pdf(here::here('plots','spot_plots_lifespan_nmf26_.pdf'), width = 20, height = 30)

###make single-variable plot with broad domain ground
ggrastr::rasterize(
plotVisium(
    spe,
    spots = TRUE,
    fill = 'V26',
    highlight = "broad.domain",
    facets = NULL,
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
        direction=1) +
    facet_wrap2(~ age, nrow = 4, ncol = 4)
)

dev.off()

# Add column in colData() for neuron
spe$neuron_cell_body <- FALSE
spe$neuron_cell_body[spe$bayesSpace_harmony_10 == "4" |
        spe$bayesSpace_harmony_10 == "7"|
        spe$bayesSpace_harmony_10 == "9"] <- TRUE

# if plotting a gene instead, add logcounts to colData(spe)
spe$ALDH1A2 <- assays(spe)$logcount["ALDH1A2",]
spe$POSTN <- assays(spe)$logcount["POSTN",]
spe$WIPF3 <- assays(spe)$logcount["WIPF3",]
spe$SFRP2 <- assays(spe)$logcount["SFRP2",]

pdf(here::here('plots', 'NMF','spot_plots_lifespan_nmf26_5_14.pdf'), width = 12, height = 22)

###make multi-variable plot (var 1=magenta, var 2=yellow, var 3=green, var 4=blue)
plotVisiumRGB(spe, pink = 'V26', yellow = 'V5', cyan = 'V14',
image=F,highlight='neuron_cell_body',values=c('gray50','black')) +
        facet_wrap2(~ age, nrow = 4, ncol = 4)

dev.off()

# Create dot plot to highlight which spatial domain nmf26 is enriched in and proportion of spots

df <-
    data.frame(spe$key, spe$sample_id, spe$bayesSpace_harmony_10)

df <- df %>%
    mutate(
        BayesSpace = case_when(
            spe.bayesSpace_harmony_10 == 1 ~ "SLM",
            spe.bayesSpace_harmony_10 == 2 ~ "ML",
            spe.bayesSpace_harmony_10 == 3 ~ "CP",
            spe.bayesSpace_harmony_10 == 4 ~ "CA3&4",
            spe.bayesSpace_harmony_10 == 5 ~ "SR",
            spe.bayesSpace_harmony_10 == 6 ~ "SGZ",
            spe.bayesSpace_harmony_10 == 7 ~ "GCL",
            spe.bayesSpace_harmony_10 == 8 ~ "SL",
            spe.bayesSpace_harmony_10 == 9 ~ "CA1",
            spe.bayesSpace_harmony_10 == 10 ~ "WM",
        )
    )

colData(spe)$BayesSpace <-
    factor(df$BayesSpace, levels = c("WM", "SLM", "SL", "ML", "CP", "SR", "SGZ",
        "CA3&4", "CA1", "GCL"))

df <- as.data.frame(colData(spe)[,c(171, 172, 95)])

# Calculate mean values of column 'V' grouped by columns 'A' and 'B'
mean_values <- aggregate(V26 ~ age_bin + BayesSpace, data = df, FUN = mean)

# Calculate the proportion of positive values in column 'V' grouped by columns 'A' and 'B'
positive_prop <- aggregate(ifelse(df$V26 > 0, 1, 0) ~ age_bin + BayesSpace, data = df, FUN = mean)

# Merge the mean and proportion dataframes
merged_df <- merge(mean_values, positive_prop, by = c("age_bin", "BayesSpace"))

# Create the dot plot
pdf(here::here('plots','dot_plots_nmf26magenta_lifespan1.pdf'))

ggplot(merged_df, aes(x = age_bin, y = BayesSpace, color = V26, size = `ifelse(df$V26 > 0, 1, 0)`)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "white", high = "magenta") +
    labs(x = "age", y = "spatial domain", color = "mean\nnmf26\nweight", size = "Proportion of spots") +
  theme_minimal()

dev.off()

# Create dot plot to highlight which spatial domain nmf5 is enriched in and proportion of spots

df2 <- as.data.frame(colData(spe)[,c(171, 172, 74)])

# Calculate mean values of column 'V' grouped by columns 'A' and 'B'
mean_values2 <- aggregate(V5 ~ age_bin + BayesSpace, data = df2, FUN = mean)

# Calculate the proportion of positive values in column 'V' grouped by columns 'A' and 'B'
positive_prop2 <- aggregate(ifelse(df2$V5 > 0, 1, 0) ~ age_bin + BayesSpace, data = df2, FUN = mean)

# Merge the mean and proportion dataframes
merged_df2 <- merge(mean_values2, positive_prop2, by = c("age_bin", "BayesSpace"))

# Create the dot plot
pdf(here::here('plots','dot_plots_nmf5yellow_lifespan1.pdf'))

ggplot(merged_df2, aes(x = age_bin, y = BayesSpace, color = V5, size = `ifelse(df2$V5 > 0, 1, 0)`)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "white", high = "yellow") +
    labs(x = "age", y = "spatial domain", color = "mean\nnmf5\nweight", size = "Proportion of spots") +
  theme_minimal()

dev.off()

# Create dot plot to highlight which spatial domain nmf14 is enriched in and proportion of spots

df3 <- as.data.frame(colData(spe)[,c(171, 172, 83)])

# Calculate mean values of column 'V' grouped by columns 'A' and 'B'
mean_values3 <- aggregate(V14 ~ age_bin + BayesSpace, data = df3, FUN = mean)

# Calculate the proportion of positive values in column 'V' grouped by columns 'A' and 'B'
positive_prop3 <- aggregate(ifelse(df3$V14 > 0, 1, 0) ~ age_bin + BayesSpace, data = df3, FUN = mean)

# Merge the mean and proportion dataframes
merged_df3 <- merge(mean_values3, positive_prop3, by = c("age_bin", "BayesSpace"))

# Create the dot plot
pdf(here::here('plots','dot_plots_nmf14green_lifespan1.pdf'))

ggplot(merged_df3, aes(x = age_bin, y = BayesSpace, color = V14, size = `ifelse(df3$V14 > 0, 1, 0)`)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "white", high = "green") +
    labs(x = "age", y = "spatial domain", color = "mean\nnmf5\nweight", size = "Proportion of spots") +
  theme_minimal()

dev.off()

# plotVisiumRGB does not produce scale legends, Make UMAP of Visium spots with color scales to add in illustrator

pdf(here::here('plots','UMAP_nmf26514_lifespan.pdf'))

ggrastr::rasterize(
plotUMAP(spe,colour_by='V26',text_by='bayesSpace_harmony_10',point_size=0.1)+
    labs(color = "nmf26\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "black", high = "magenta") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

ggrastr::rasterize(
plotUMAP(spe,colour_by='V5',text_by='bayesSpace_harmony_10',point_size=0.1)+
    labs(color = "nmf5\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "black", high = "yellow") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

ggrastr::rasterize(
plotUMAP(spe,colour_by='V14',text_by='bayesSpace_harmony_10',point_size=0.1)+
    labs(color = "nmf5\nweight", max.overlaps = Inf) +
        scale_color_gradient(low = "black", high = "green") +
        theme(axis.ticks=element_blank(),axis.text=element_blank())
)

dev.off()
