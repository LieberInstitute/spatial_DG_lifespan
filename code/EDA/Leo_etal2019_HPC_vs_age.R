####################################################################
# Checking Jaffe et al, 2020 Nat. Neuro. for age changes
####################################################################

##### libraries
library(SummarizedExperiment)
library(jaffelab)
library(sva)
library(limma)
library(edgeR)
library(recount)
library(RColorBrewer)
library(ggplot2)

######################
## load data #########
######################

load("rse_gene_unfiltered.Rdata")

rse_HPC_gene <- rse_gene[, which(rse_gene$Region == "HIPPO")]
rse_HPC_gene <- rse_HPC_gene[, which(rse_HPC_gene$Dx == "Control")]

rownames(rse_HPC_gene) <- rowData(rse_HPC_gene)$Symbol

## plots
exampleIndex = c("PROX1", "CALB1", "GAD1", "LAMP5", "SLC17A7", "CD74", "CCK", "GAD2", "POSTN", "DCX")

# Extract RPKM values
rpkm_values <- assays(rse_HPC_gene)$rpkm

# Extract age information
age <- rse_HPC_gene$Age

# Loop through each gene and create a plot
# Create a list to store the individual plots
gene_plots <- list()

for (i in exampleIndex) {
  # Extract RPKM values for the current gene
  gene_rpkm_values <- rpkm_values[i, ]

  # Create a data frame for the current gene
  gene_data <- data.frame(Gene = i, RPKM = gene_rpkm_values, Age = age)

  # Create a plot for the current gene
  gene_plot <- ggplot(gene_data, aes(x = Age, y = RPKM)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    labs(x = "Age", y = "RPKM Values") +
    ggtitle(paste("RPKM Values vs. Age for", i)) +
      theme(
        axis.line = element_line(colour = "black"),
        text = element_text(colour = 'black', size = 10),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
    theme_classic()

  # Store the plot in the list
  gene_plots[[i]] <- gene_plot
}

# Print or display the individual plots
pdf("ageEffect_HPC_example_genes.pdf")
for (i in exampleIndex) {
  print(gene_plots[[i]])
}

dev.off()

