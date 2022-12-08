####################################
# spatial_DG_lifespan project
# trucate sce_sestan object clusters
# Anthony Ramnauth, Nov 05 2022
####################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(dplyr)
    library(sessioninfo)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "QCed_sce_sestan_DG.rds"))

unique(sce$cluster)

df <- as.data.frame(colData(sce))

df <- df %>%
  mutate(Cell_Type = case_when(
    grepl("CA3", cluster) ~ "CA3_N",
    grepl("DG MC", cluster) ~ "Mossy",
    grepl("CA1 d", cluster) ~ "CA1_d_N",
    grepl("CA1 v", cluster) ~ "CA1_v_N",
    grepl("CA2", cluster) ~ "CA2_N",
    grepl("SGCZ", cluster) ~ "GC_1",
    grepl("PDLIM5", cluster) ~ "GC_2",
    grepl("InN LAMP5", cluster) ~ "InN_LAMP5",
    grepl("InN VIP", cluster) ~ "InN_VIP",
    grepl("OTOF", cluster) ~ "InN_SST",
    grepl("ADAMTS12", cluster) ~ "InN_SST",
    grepl("EPB41L4A", cluster) ~ "InN_SST",
    grepl("NPY", cluster) ~ "InN_NPY",
    grepl("InN PVALB", cluster) ~ "InN_PV",
    grepl("InN NR2F2", cluster) ~ "InN_NR2F2",
    grepl("InN LHX6", cluster) ~ "InN_LHX6",
    grepl("InN MEIS2", cluster) ~ "InN_MEIS2",
    grepl("CR RELN NDNF", cluster) ~ "Cajal_Ret",
    grepl("VLMC", cluster) ~ "Vasc_LM",
    grepl("GRIA4", cluster) ~ "OPC_1",
    grepl("EGR1", cluster) ~ "OPC_2",
    grepl("aSMC", cluster) ~ "Artl_S_Muscle",
    grepl("CPXM2", cluster) ~ "Oligo_2",
    grepl("OPALIN", cluster) ~ "Oligo_1",
    grepl("Micro", cluster) ~ "Microglia",
    grepl("PC CLDN5", cluster) ~ "Pericyte",
    grepl("Endo", cluster) ~ "Endoth",
    grepl("vSMC", cluster) ~ "Vasc_S_Muscle",
    grepl("Macro", cluster) ~ "Macrophage",
    grepl("aEndo", cluster) ~ "Endo",
    grepl("COP GPR17", cluster) ~ "COP",
    grepl("T SKAP1", cluster) ~ "T_cell",
    grepl("GFAP", cluster) ~ "Astro_1",
    grepl("CHRDL1", cluster) ~ "Astro_2",
    grepl("Myeloid", cluster) ~ "Myeloid",
  ))

colData(sce)$Cell_Type <- df$Cell_Type

saveRDS(sce, file = here::here("processed-data","sce", "sce_sestan_DG_final.rds"))
