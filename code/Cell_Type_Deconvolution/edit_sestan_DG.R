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
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "QCed_sce_sestan_DG.rds"))

unique(sce$cluster)

df <- as.data.frame(colData(sce))

df <- df %>%
  mutate(Cell_Type = case_when(
    grepl("CA3", cluster) ~ "CA3_N",
    grepl("DG MC", cluster) ~ "Mossy",
    grepl("CA1", cluster) ~ "CA1_N",
    grepl("CA2", cluster) ~ "CA2_N",
    grepl("PROX1", cluster) ~ "GC",
    grepl("InN LAMP5", cluster) ~ "InN_LAMP5",
    grepl("InN VIP", cluster) ~ "InN_VIP",
    grepl("InN SST", cluster) ~ "InN_SST",
    grepl("InN PVALB", cluster) ~ "InN_PV",
    grepl("InN NR2F2", cluster) ~ "InN_NR2F2",
    grepl("InN LHX6", cluster) ~ "InN_LHX6",
    grepl("InN MEIS2", cluster) ~ "InN_MEIS2",
    grepl("CR RELN NDNF", cluster) ~ "Cajal_Retz",
    grepl("VLMC", cluster) ~ "VLMC",
    grepl("PDGFRA", cluster) ~ "OPC",
    grepl("aSMC", cluster) ~ "SMC",
    grepl("CPXM2", cluster) ~ "Oligo",
    grepl("OPALIN", cluster) ~ "Oligo",
    grepl("Micro", cluster) ~ "Microglia",
    grepl("PC CLDN5", cluster) ~ "Pericyte",
    grepl("Endo", cluster) ~ "Endoth",
    grepl("vSMC", cluster) ~ "SMC",
    grepl("Macro", cluster) ~ "Macro",
    grepl("aEndo", cluster) ~ "Endoth",
    grepl("COP GPR17", cluster) ~ "COP",
    grepl("T SKAP1", cluster) ~ "T_Cell",
    grepl("GFAP", cluster) ~ "Astro_1",
    grepl("CHRDL1", cluster) ~ "Astro_2",
    grepl("Myeloid", cluster) ~ "Myeloid",
  ))

  #as.data.frame(table(df$Cell_Type))
  #         Var1  Freq
  # 1    Astro_1  8607
  # 2    Astro_2  4331
  # 3      CA1_N    33
  # 4      CA2_N    25
  # 5      CA3_N  1403
  # 6        COP    65
  # 7     Endoth  1111
  # 8         GC 24613
  # 9  InN_LAMP5   789
  # 10  InN_LHX6    65
  # 11 InN_MEIS2    25
  # 12 InN_NR2F2   442
  # 13    InN_PV   255
  # 14   InN_SST   469
  # 15   InN_VIP   990
  # 16     Macro   108
  # 17 Microglia  6340
  # 18     Mossy   280
  # 19   Myeloid    15
  # 20     Oligo 41277
  # 21       OPC 10558
  # 22  Pericyte   145
  # 23       SMC   362
  # 24    T_Cell   147
  # 25      VLMC   298
  
colData(sce)$Cell_Type <- df$Cell_Type


saveRDS(sce, file = here::here("processed-data","sce", "sce_sestan_DG_final.rds"))

