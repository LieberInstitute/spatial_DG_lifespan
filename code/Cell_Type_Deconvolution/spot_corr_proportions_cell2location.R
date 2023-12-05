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
