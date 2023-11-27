########################################################
# spatial_DG_lifespan project
# Approx. spot correlation of cellular neighborhood info
# Anthony Ramnauth, Nov. 13 2023
########################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggcorrplot)
    library(ggsignif)
    library(ggh4x)
    library(crawdad)
    library(Seurat)
    library(sessioninfo)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

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

spe <- spe[, order(spe$age)]

# Convert mean abundances to proportions

convert_to_proportions <- function(row) {
  proportions <- row / sum(row)
  return(proportions)
}


spe_props <- as.matrix(colData(spe)[, c(44:68)])
# Apply the function to each row of the data frame
spe_props <- apply(spe_props, 1, convert_to_proportions)
# Convert the result back to a data frame
spe_props <- as.data.frame(t(spe_props))
for (col in 1:ncol(spe_props)){
    colnames(spe_props)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(spe_props)[col])
}

# Compute correlation matrix
cor_props = cor(spe_props)

colors <- c("#91a28c","white","#8f2c37")

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "cell2loc_spot_corr.pdf"), width = 9.5, height = 9.5)

ggcorrplot(cor_props, hc.order = F,
   outline.color = "white",
   tl.srt = 60,
   tl.cex = 18,
   lab_size = 7,
   colors = colors)+
 theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
     axis.text = element_text(size = 12),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 16,face="bold"),
    legend.text=element_text(size = 16),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm')) +
 coord_fixed()+
 ggtitle("Cell type proportion spot correlation")+theme(plot.title = element_text(size=22,face="bold"))

dev.off()

# Use CRAWDAD package to compute neighborhood of dominant cell types for spots.

# Find dominant celltype for each spot
spe_dom <- spe_props %>%
 rowwise() %>%
 mutate(row_max = names(.)[which.max(c_across(everything()))])

colData(spe)$dominant_cell_types <- spe_dom$row_max

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
colData(spe)$x <-
    colData(spe)$array_row + auto_offset_row[spe$sample_id]
colData(spe)$y <- colData(spe)$array_col

colnames(spe) <- spe$key

# Convert to Seurat format for CRAWDAD

DG_seu <- CreateSeuratObject(
      counts=counts(spe),
      meta.data=data.frame(colData(spe)),
      project="DG_lifespan")

craw_obj <- crawdad:::seuratToSF(DG_seu,
   coms  = "dominant_cell_types",
    posIDs = c("x", "y"),
    verbose = TRUE
    )

scales <- seq(100, 1000, by=100)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(craw_obj,
                          scales = scales,
                          perms = 3,
                          ncores = 7,
                          seed = 1,
                          verbose = FALSE)

## find trends, passing background as parameter
results <- crawdad::findTrends(craw_obj,
                        dist = 100,
                        shuffle.list = shuffle.list,
                        ncores = 7,
                        verbose = FALSE,
                        returnMeans = FALSE)

## convert results to data.frame
dat <- crawdad::meltResultsList(results, withPerms = T)

# Save dataframe

saveRDS(dat, file = here::here("processed-data", "Cell_Type_Deconvolution", "crawdad_results.rds"))

## multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

pdf(file = here::here("plots", "Cell_Type_Deconvolution", "crawdad_z_scale.pdf"))

vizColocDotplot(dat, reorder = TRUE, zsig.thresh = zsig, zscore.limit = zsig*2) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

dev.off()

