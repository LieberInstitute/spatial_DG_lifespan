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


## load data (can be downloaded here: https://eqtl.brainseq.org/phase2/ )

load("rse_gene_unfiltered.Rdata")

rse_HPC_gene <- rse_gene[, which(rse_gene$Region == "HIPPO")]
rse_HPC_gene <- rse_HPC_gene[, which(rse_HPC_gene$Dx == "Control")]

rownames(rse_HPC_gene) <- rowData(rse_HPC_gene)$Symbol

## plots
exampleIndex = c("PROX1", "CALB1", "GAD1", "LAMP5", "SLC17A7", "CD74", "CCK", "GAD2", "POSTN", "DCX")

x <- edgeR::cpm(edgeR::calcNormFactors(rse_HPC_gene), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(rse_HPC_gene)))

# Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(rse_HPC_gene)

# Store the log normalized counts on the RangedSummarizedExperiment object
assays(rse_HPC_gene)$logcounts <- x

dim(rse_HPC_gene)

# Extract age information
age <- rse_HPC_gene$Age

# Loop through each gene and create a plot
# Create a list to store the individual plots
gene_plots <- list()

for (i in exampleIndex) {
  # Extract RPKM values for the current gene
  gene_values <- x[i, ]

  # Create a data frame for the current gene
  gene_data <- data.frame(Gene = i, normexpr = gene_values, Age = age)

  # Create a plot for the current gene
  gene_plot <- ggplot(gene_data, aes(x = Age, y = normexpr)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    ylim(0,NA) +
    labs(x = "Age", y = "Normalized Expression") +
    ggtitle(paste("Normalized Expression vs. Age for", i)) +
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

