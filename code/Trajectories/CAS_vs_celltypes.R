################################################
# spatial_DG_lifespan project
# Checking if cell type composition predicts CAS
# Anthony Ramnauth, Sept 20 2023
################################################

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(dplyr)
    library(stats)
    library(ggplot2)
})

# Load SPE with cell type mean abundances
spe_cell <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

# order spe observations according to age
spe_cell <- spe_cell[, order(spe_cell$age)]

# Make data frame with  cell abundances and use key as common rownames
colnames(spe_cell) <- spe_cell$key

cell_abundances <- as.data.frame(colData(spe_cell)[, 44:68])

# Get dominant cell types per Visium spot for plots

cell_props <- as.data.frame(colData(spe_cell)[, c(44:68)])
for (col in 1:ncol(cell_props)){
    colnames(cell_props)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(cell_props)[col])
}

# Find dominant celltype for each spot
p <- cell_props %>%
 rowwise() %>%
 mutate(row_max = names(.)[which.max(c_across(everything()))])

colData(spe_cell)$dominant_cell_types <- p$row_max

# Set colors
cell_colors <- c("Oligo" = "plum4","Microglia" = "tan3", "Macro" = "tan4","OPC" = "goldenrod",
    "InN_LAMP5" = "springgreen", "InN_VIP" = "green1", "InN_SST" = "springgreen2", "InN_PV" = "green",
    "InN_NR2F2" = "green2", "InN_LHX6" = "springgreen", "InN_MEIS2" = "springgreen3",
    "VLMC" = "firebrick1", "Pericyte" = "red2", "Endoth" = "red", "SMC" = "firebrick",
    "T_Cell" = "brown1", "Myeloid" = "tan", "COP" = "goldenrod4", "GC" = "blue2","CA3_N" = "navy",
    "Mossy" = "blue1", "CA1_N" = "blue3","CA2_N" = "blue4", "Astro_1" = "magenta", "Astro_2" = "yellow3")

# Load SPE with cell type mean abundances
spe_age <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_CAS_BayesSpace_spe.rds"))

# Make data frame with CAS and use key as common rownames
colnames(spe_age) <- spe_age$key

CAS <- as.data.frame(as.data.frame(colData(spe_age)[, 44:45]))

# Step 1: Calculate the 1st principal component of 'sup'
pca_result <- prcomp(cell_abundances, center = TRUE, scale = TRUE)

jaffelab::getPcaVars(pca_result)[seq_len(25)]
#[1] 14.500 12.300  8.130  7.090  6.420  5.220  3.930  3.700  3.640  3.470  3.160  3.000  2.940  2.810  2.640  2.480  2.390  2.220
#[19]  2.040  1.960  1.700  1.570  1.200  0.759  0.719

# Step 2: Create a new data frame with 'blah' and the 1st principal component
new_data <- data.frame(CAS = CAS$CAS,
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3],
    PC4 = pca_result$x[, 4]
)

# Fit a linear model to the log-transformed data
fit_cas_abund_PC1 <- lm(CAS ~ PC1, data = new_data)
fit_cas_abund_PC2 <- lm(CAS ~ PC2, data = new_data)
fit_cas_abund_PC3 <- lm(CAS ~ PC3, data = new_data)
fit_cas_abund_PC4 <- lm(CAS ~ PC4, data = new_data)


# Extract the coefficients and R2 value
intercept_abund_PC1 <- coef(fit_cas_abund_PC1)[1]
slope_abund_PC1 <- coef(fit_cas_abund_PC1)[2]
r2_abund_PC1 <- summary(fit_cas_abund_PC1)$r.squared

intercept_abund_PC2 <- coef(fit_cas_abund_PC2)[1]
slope_abund_PC2 <- coef(fit_cas_abund_PC2)[2]
r2_abund_PC2 <- summary(fit_cas_abund_PC2)$r.squared

intercept_abund_PC3 <- coef(fit_cas_abund_PC3)[1]
slope_abund_PC3 <- coef(fit_cas_abund_PC3)[2]
r2_abund_PC3 <- summary(fit_cas_abund_PC3)$r.squared

intercept_abund_PC4 <- coef(fit_cas_abund_PC4)[1]
slope_abund_PC4 <- coef(fit_cas_abund_PC4)[2]
r2_abund_PC4 <- summary(fit_cas_abund_PC4)$r.squared

# plot CAS vs PCs

pdf(file = here::here("plots", "Trajectories", "celltype_abundances_PC_vs_CAS.pdf"), width = 9, height = 5)

ggplot(new_data,
    aes(
        x = PC1,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("1st Principal Component of mean cell-type abundances") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = -23, y = 1,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_abund_PC1, r2_abund_PC1),
        size = 4) +
    annotate("text", x = -25, y = .85,
        label = sprintf("R2 = %.2f", r2_abund_PC1), size = 4)

ggplot(new_data,
    aes(
        x = PC2,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("2nd Principal Component of mean cell-type abundances") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 15, y = 1,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_abund_PC2, r2_abund_PC2),
        size = 4) +
    annotate("text", x = 15, y = .85,
        label = sprintf("R2 = %.2f", r2_abund_PC2), size = 4)

ggplot(new_data,
    aes(
        x = PC3,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("3rd Principal Component of mean cell-type abundances") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 10, y = 1,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_abund_PC3, r2_abund_PC3),
        size = 4) +
    annotate("text", x = 10, y = .85,
        label = sprintf("R2 = %.2f", r2_abund_PC3), size = 4)

ggplot(new_data,
    aes(
        x = PC4,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("4th Principal Component of mean cell-type abundances") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = -12, y = 1,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_abund_PC4, r2_abund_PC4),
        size = 4) +
    annotate("text", x = -12, y = .85,
        label = sprintf("R2 = %.2f", r2_abund_PC4), size = 4)

dev.off()

# Convert to cell type proportions per spot and rerun

# Convert mean abundances to proportions

convert_to_proportions <- function(row) {
  proportions <- row / sum(row)
  return(proportions)
}

# Apply the function to each row of the data frame
cell_proportions <- apply(cell_abundances, 1, convert_to_proportions)
# Convert the result back to a data frame
cell_proportions <- as.data.frame(t(cell_proportions))

# Step 1: Calculate the 1st principal component of 'sup'
pca_result2 <- prcomp(cell_proportions, center = TRUE, scale = TRUE)

jaffelab::getPcaVars(pca_result2)[seq_len(25)]
#[1] 1.34e+01 8.66e+00 8.45e+00 7.21e+00 6.22e+00 4.73e+00 4.15e+00 3.99e+00 3.93e+00 3.70e+00 3.64e+00 3.45e+00 3.24e+00 3.13e+00
#[15] 3.06e+00 2.90e+00 2.69e+00 2.60e+00 2.42e+00 2.21e+00 2.02e+00 1.74e+00 1.46e+00 1.01e+00 1.33e-28

# Step 2: Create a new data frame with 'blah' and the 1st principal component
new_data2 <- data.frame(CAS = CAS$CAS,
    PC1 = pca_result2$x[, 1],
    PC2 = pca_result2$x[, 2],
    PC3 = pca_result2$x[, 3],
    PC4 = pca_result2$x[, 4]
)

# Fit a linear model to the log-transformed data
fit_cas_prop_PC1 <- lm(CAS ~ PC1, data = new_data2)
fit_cas_prop_PC2 <- lm(CAS ~ PC2, data = new_data2)
fit_cas_prop_PC3 <- lm(CAS ~ PC3, data = new_data2)
fit_cas_prop_PC4 <- lm(CAS ~ PC4, data = new_data2)

# Extract the coefficients and R2 value
intercept_prop_PC1 <- coef(fit_cas_prop_PC1)[1]
slope_prop_PC1 <- coef(fit_cas_prop_PC1)[2]
r2_prop_PC1 <- summary(fit_cas_prop_PC1)$r.squared

intercept_prop_PC2 <- coef(fit_cas_prop_PC2)[1]
slope_prop_PC2 <- coef(fit_cas_prop_PC2)[2]
r2_prop_PC2 <- summary(fit_cas_prop_PC2)$r.squared

intercept_prop_PC3 <- coef(fit_cas_prop_PC3)[1]
slope_prop_PC3 <- coef(fit_cas_prop_PC3)[2]
r2_prop_PC3 <- summary(fit_cas_prop_PC3)$r.squared

intercept_prop_PC4 <- coef(fit_cas_prop_PC4)[1]
slope_prop_PC4 <- coef(fit_cas_prop_PC4)[2]
r2_prop_PC4 <- summary(fit_cas_prop_PC4)$r.squared

# plot CAS vs PCs

pdf(file = here::here("plots", "Trajectories", "celltype_proportions_PC_vs_CAS.pdf"), width = 9, height = 5)

ggplot(new_data2,
    aes(
        x = PC1,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("1st Principal Component of cell-type proportions (13.4%)") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = -8, y = 2,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_prop_PC1, r2_prop_PC1),
        size = 4) +
    annotate("text", x = -8, y = 1.85,
        label = sprintf("R2 = %.2f", r2_prop_PC1), size = 4)

ggplot(new_data2,
    aes(
        x = PC2,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("2nd Principal Component of cell-type proportions (8.66%)") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 10, y = 2,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_prop_PC2, r2_prop_PC2),
        size = 4) +
    annotate("text", x = 10, y = 1.85,
        label = sprintf("R2 = %.2f", r2_prop_PC2), size = 4)

ggplot(new_data2,
    aes(
        x = PC3,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("3rd Principal Component of cell-type proportions (8.45%)") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 5, y = 2,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_prop_PC3, r2_prop_PC3),
        size = 4) +
    annotate("text", x = 5, y = 1.85,
        label = sprintf("R2 = %.2f", r2_prop_PC3), size = 4)

ggplot(new_data2,
    aes(
        x = PC4,
        y = CAS,
        color = factor(spe_cell$dominant_cell_types)
        )) +
    geom_point(size = 0.5, alpha = 0.3) +
    scale_color_manual(values = cell_colors) +
    xlab("4th Principal Component of cell-type proportions (7.21%)") +
    ylab("CAS") +
    labs(color = "Dominant cell type") +
    theme_bw() +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black",linetype = 'dashed',size = 0.5) +
    theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    annotate("text", x = 9, y = 2,
        label = sprintf("Equation: y = %.2f + %.2fx", intercept_prop_PC4, r2_prop_PC4),
        size = 4) +
    annotate("text", x = 9, y = 1.85,
        label = sprintf("R2 = %.2f", r2_prop_PC4), size = 4)

dev.off()

