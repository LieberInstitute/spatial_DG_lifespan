####################################
# spatial_DG_lifespan project
# trucate sce_sestan object clusters
# Anthony Ramnauth, Nov 05 2022
####################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(here)
    library(dplyr)
    library(sessioninfo)
})

# load saved sce object

sce <- readRDS(here::here("processed-data", "sce", "sce_sestan_DG.rds"))

df <- as.data.frame(colData(sce))

df <- df %>%
  mutate(Cell_Type = case_when(
    grepl("CA3", cluster) ~ "CA3_N",
    grepl("EC L", cluster) ~ "EC_N",
    grepl("DG MC", cluster) ~ "Mossy",
    grepl("CA1", cluster) ~ "CA1_N",
    grepl("SUB distal", cluster) ~ "SUB_N",
    grepl("SUB proximal", cluster) ~ "SUB_N",
    grepl("SGCZ", cluster) ~ "GC_1",
    grepl("PDLIM5", cluster) ~ "GC_2",
    grepl("InN LAMP5", cluster) ~ "InN_LAMP5",
    grepl("InN VIP", cluster) ~ "InN_VIP",
    grepl("InN SST", cluster) ~ "InN_SST",
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
    grepl("Macro", cluster) ~ "Macrocyte",
    grepl("aEndo", cluster) ~ "Endo",
    grepl("COP GPR17", cluster) ~ "COP",
    grepl("T SKAP1", cluster) ~ "T_cell",
    grepl("GFAP", cluster) ~ "Astro_1",
    grepl("CHRDL1", cluster) ~ "Astro_2",
    grepl("Myeloid", cluster) ~ "Myeloid",
  ))

colData(sce)$Cell_Type <- df$Cell_Type

saveRDS(sce, file = here::here("sce_objects", "sce_sestan_DG.rds"))
