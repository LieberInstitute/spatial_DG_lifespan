Layer-level documentation
===========================================

## Common options

* `Model results`: the statistical modeling results to use. We computed three different types of models:
  1. `enrichment`: one spatial domain against all the other domains. Results in t-statistics. `GCL`, for example, represents enriched and depleted genes in the `GCL` (`next_Ab`) spatial domain against all the other domains based on this enrichment modeling method.
  2. `pairwise`: one spatial domain against another one. Results in t-statistics with two-sided p-values. GAD2 in `SGZ` vs `GCL`, for example, compares gene expression of GAD between the `SGZ` (`Abeta`) and `GCL` spatial domain. 
  3. `anova`: changes among the spatial domains (adjusting for the mean expression) using the data from all 9 spatial domains (Choroid Plexus is dropped from pseudobulked data).

## Reduced dim

In this panel you can visualize the layer-level data (`spe_pseudo`) across reduced dimensionality representations derived from the gene expression data from the layer-level pseudo-bulked data. Select which dimensionality reduction method (here we provide PCA). Then use `Color by` to choose which variable to color data by, which can be useful to identify groups of pseudo-bulked samples. The options include e.g.,

* `age`: age of death of n = 16 donors.
* `BayesSpace`: Spatial domain of n = 16 donors.
* `Diagnosis`: clinical diagnosis of n = 16 donors.
* `ncells`: number of spots that were combined when pseudo-bulking.
* `pmi`: post-mortem interval of n = 3 donors with AD.
* `race`: race of n = 3 donors with AD.
* `rin`: RNA integrity number of n = 3 donors with AD.
* `sample_id`: sample identifier of n = 16 samples from n = 16 donors.
* `sex`: sex of n = 16 donors.


```{r}
## Reproduce locally with
scater::plotReducedDim(sce_pseudo)
```

## Model boxplots

This tab allows you to make a boxplot of the `logcounts` gene expression from the spatial domain-level data (`spe_pseudo`) for a given `gene`; you can search your gene by typing either the symbol or the Ensembl gene ID. The model result information displayed in the title of the plot is based on which `model results` you selected and whether you are using the short title version or not (controlled by a checkbox). We provide two different color scales you can use: the color blind friendly `viridis` as well as a custom one we used for the `paper`. Through the `Model test` selector, you can choose which particular comparison to display.


Below the plot you can find the subset of the table of results  (`sig_genes`), sort the table by the different columns, and download it as a CSV if you want. For more details about what each of these columns mean, check the [`spatialLIBD` vignette documentation](http://LieberInstitute.github.io/spatialLIBD/articles/spatialLIBD.html#extract-significant-genes).

```{r}
## Reproduce locally with
spatialLIBD::layer_boxplot()
```

## Gene Set Enrichment

This tab allows you to upload a CSV file that has a particular format as illustrated [in this example file](https://github.com/LieberInstitute/spatialLIBD/blob/master/data-raw/asd_sfari_geneList.csv). This CSV file should contain:

* one column per gene set of interest labeled as column names on the first row,
* no row names, 
* and human Ensembl gene IDs as values in the cells. 

Once you have uploaded a CSV file following this specific format, you can then check if the genes on each of your gene sets are enriched among the statistics from `model results` (`enrichment`, etc) that have a false discovery rate (FDR) adjusted p-value less than `FDR cutoff` (0.1 by default).

Similar to the `Model boxplots` tab, you can interact with the results table or download it.

```{r}
## Reproduce locally with
spatialLIBD::gene_set_enrichment()
spatialLIBD::gene_set_enrichment_plot()
```

## Spatial registration

If you have a single nucleus or single cell RNA-sequencing (snRNA-seq)  (scRNA-seq) dataset, you might group your cells into clusters. Once you do, you could compress the data by pseudo-bulking (like we did to go from `spe` to `sce_pseudo`). You could then compute `enrichment` (`pairwise`, `anova`) statistics for your cell clusters. If you do so, you can then upload a specially formatted CSV file just like the one in [this example file](https://github.com/LieberInstitute/spatialLIBD/blob/master/data-raw/tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer.csv). This file has:

* column names,
* human Ensembl gene IDs as the row names (first column, no name for the column),
* statistics (numeric values) for the cells.

Once you have uploaded a CSV file following this specific format, you can then assess whether the correlation between your statistics and the ones from our spatial domains for the subset of genes (Ensembl ids) present in both. The resulting heatmap and interactive correlation matrix (which again you can interact with and download) can be useful if you are in the process of labeling your sn/scRNA-seq clusters or simply want to compare them against the spatial domain-specific data we have provided. This can also be used for new spatially-resolved transcriptomics datasets.

Finally, you can change the `Maximum correlation` for visualization purposes on the heatmap as it will change the dynamic range for the colors.

```{r}
## Reproduce locally with
spatialLIBD::layer_stat_cor()
spatialLIBD::layer_stat_cor_plot()
```
