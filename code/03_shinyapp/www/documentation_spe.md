Spot-level documentation
======================================

## Common options

* `Samples to plot`: which sample to plot on the tabs that do not have _grid_ on their name. Samples are named by donor-- in practice, tissue was split across several Visium slides and computationally stitched together to include the entire nucleus accumbens for each donor. [TODO: describe more about donors]
* `Image name`: the name of the background image to use. You can edit this image on the `Edit image` tab.
  - `edited_imaged`: shows your edited image from the Edit image tab.
  - `lowres`: shows low resolution images of raw immunofluorescence data before image processing such as spectral unmixing and object segmentation.
* `Discrete variable to plot`: which discrete variable (typically with the cluster labels) to visualize. We include the clusters:
  - `BayesSpace`: Unsupervised spatial clustering of the spots.
  - `ManualAnnotation`: your own manual annotation of the spots.
* `Continuous variable to plot`: which gene or continuous variable (such as the cell count, the ratio of the mitochondrial chromosome expression) to visualize in the gene tabs as well as on the `clusters (interactive)` tab. Details:
  - `sum_umi`: sum of UMI counts across a spot.
  - `sum_gene`: number of genes with non-zero counts in a spot.
  - `expr_chrM`: sum of chrM counts in a spot.
  - `expr_chrM_ratio`: ratio of `expr_chrM / sum_umi`.
* `Spot transparency level`: the transparency of the spots in the visualizations. It can be useful if the spot colors are blocking the background image.
* `Spot point size`: the size of the spots in the visualizations.
* `Minimum count value`: Values from the selected `continuous variable to plot` at or below this threshold will not be displayed.
* `Gene color scale`: Whether to use the color blind friendly palette (`viridis`) or to use a custom palette that we used for our `paper`. Other options from the [viridisLite R package](https://sjmgarnier.github.io/viridisLite/reference/viridis.html#details) are also supported.
* `Gene color direction`: whether colors should be ordered from darkest to lightest or in the reverse direction.

We will cover the download button and upload CSV options at the end of this document.

## Clusters (static)

Displays the selected cluster variable (from `discrete variable to plot`) for the given sample (from `samples to plot`). You can also choose to see the background image and clusters side to side.

```{r}
## Reproduce locally with
spatialLIBD::vis_clus()
```

## Clusters (interactive)

Displays a 1,200 by 1,200 pixels interactive plot area with a matrix of 2 by 2 plots. The top row shows the data at the spot-level with the histology information in the background. The bottom row shows the spot-level data at a reduced dimension space (PCA, TSNE, UMAP, etc). The left column shows the selected gene or continuous variable, while the right column shows the selected cluster or discrete variable. The four plots are linked to each other such that if you use the lasso selector (mouse over to the top right of the interactive area to select it) in a single plot, the other 3 will get updated to highlight the same selection of points.

This panel allows you to look at the results from a given clustering approach and combine that information with the expression of a given gene to visualize the spot-level data both in the spatial resolution as well as a reduced dimensionality space from the expression of the most variable genes.

Once you have selected spots of interest, at the bottom of the tab there is a text box where you can enter your manual annotations. This overwrites `spe$ManualAnnotation` which is why you need to confirm doing so by clicking the button `Label selected points (from lasso) with manual annotation`. You can then change the `clusters to plot` option to `ManualAnnotation` to see your new spot labels.

Given the amount of data being displayed, this tab consumes quite a bit of resources.

Note that there can be a maximum of 36 unique manual annotations (unique values in `spe$ManualAnnotation`) before we run out of colors to visualize them.

## Clusters grid (static)

This tab is similar to `clusters (static)` with the difference being that you can visualize a subset or even all samples at once. Select the samples you want to visualize, then specify the grid on which they will be plotted (the number of rows and columns), and click the `upgrade grid plot` button. We ask that you use this button every time you want to update the plot as it can take a few seconds for this update to complete.

This particular tab is useful if you have computationally defined some clusters using the data of all your samples and you want to visualize the resulting clusters across multiple (or even all the) samples at the same time.

```{r}
## Reproduce locally with
spatialLIBD::vis_grid_clus()
```

## Gene (static)

This tab is similar to `clusters (static)` but instead of displaying discrete values (like clusters), it displays continuous values such as the gene expression of a given gene or the number of cells per spot. By default, spots whose value is below or at 0 are not shown, which makes it easier for you to distinguish points with low values from those below the threshold of your interest (controlled by `minimum count value`). The points can be colored in two different color scales. You can also choose to see the background image and continuous values side to side.

```{r}
## Reproduce locally with
spatialLIBD::vis_gene()
```

## Gene (interactive)

This tab shows a single interactive plot. It is similar to `clusters (interactive)` as in you can use the lasso selector (mouse over top right to find it) to select spots and them label them using the text box at the bottom with the corresponding button. Unlike `clusters (interactive)`, this version includes checkboxes at the top for each of the unique values of the selected `discrete variable to plot` such that you can subset the spots to those present on a given cluster. 

Note that if you have `discrete variable to plot` toggled to the `ManualAnnotation` option and update the spot-level information with the text box and button at the bottom of the tab, then it will re-load the interactive visualization and you will lose your selection of points.

## Gene grid (static)

This is the equivalent of `clusters grid (static)` but for continuous variables, just like `gene (static)` was the equivalent of `clusters (static)`.

```{r}
## Reproduce locally with
spatialLIBD::vis_grid_gene()
```

## Edit image

This panel shows you all the options we have for manipulating the colors and properties of the selected background image. These results will be available under the `edited_image` option in the image chooser menu. You can perform sequential manipulations of the image, though it will get hard to keep track of all the changes you made.

```{r}
## Reproduce locally with
spatialLIBD::img_update_all()
```

## Saving and uploading your `spe$ManualAnnotation` results

Beyond visualizing the data, the main goal of this section of `spatialLIBD` is to enable you to label spots. That is done only through the `interactive` tabs as previously described. Your manual annotations are saved under `spe$ManualAnnotation` and you can save them for future use using the main menu download button. If you click this button, it will prompt you to save a CSV file. We recommend that you keep selected the checkbox `Drop NA layer entries in the CSV file?` which means that your CSV file will be smaller and will potentially lead to less conflicts with other CSV files you make. That is, you will likely avoid re-labeling the same spot with different values in two or more of these CSV files. This is particularly useful if you plan to work on one sample at a time, save your results, and then merge them all into a single CSV file.

These CSV files with your manual annotations can be re-uploaded to `spatialLIBD`. You can notice this by choosing `ManualAnnotation` under `clusters to plot` and using any of the `clusters` tabs. If you upload more than one CSV, any values you have under `spe$ManualAnnotation` will be overwritten if the spot is present in your CSV file. Thus, if you followed the recommended workflow of saving one CSV file per sample, you can then upload them all sequentially and merge them together into a single CSV file to simplify your work later.

In summary, the order in which you re-upload the CSV files matters as newer uploads will overwrite any duplicated spots from previous CSV files.

We also recommend saving your work often in case you lose connection to `spatialLIBD`. Though you could always run this website locally by using the following command:

```{r}
## Run this web application locally with:
spatialLIBD::run_app()

## You will have more control about the length of the session and memory usage.
## See http://research.libd.org/spatialLIBD/reference/run_app.html#examples.
## Also see https://github.com/LieberInstitute/spatial_NAc/tree/main/code/06_deploy_app.

## You could also use spatialLIBD::run_app() to visualize your
## own data given some requirements described
## in detail in the package vignette documentation
## at http://research.libd.org/spatialLIBD/.
```

This will require about 1.8GB of RAM to run on the server side, though potentially more, specially when using the `clusters (interactive)` tab.
