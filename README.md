# spatial_DG_lifespan
Spatial transcriptomics (Visium) in dentate gyrus of hippocampus over the lifespan

spatial_DG_lifespan
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/377886452.svg)](https://zenodo.org/badge/latestdoi/377886452)

## Overview

<img src="https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/img/study_overview.png" align="left" width="300px" />

Welcome to the Spatial_DG_lifespan project! It is composed of:
* a shiny web application that we are hosting at [https://libd.shinyapps.io/Lifespan_DG/](https://libd.shinyapps.io/Lifespan_DG/) that can handle a limited set of concurrent users,  
* and a research article with the scientific knowledge we drew from this dataset. The analysis code for our project is available [here](https://github.com/LieberInstitute/spatial_DG_lifespan) and the high quality figures for the manuscript are available through Figshare. 

The web application allows you to browse the LIBD human lifespan dentate gyrus (DG) spatial transcriptomics data generated with the 10x Genomics Visium platform. Please check the manuscript or bioRxiv pre-print for more details about this research.
If you tweet about this website, the data or the R package please use the #spatialLIBD hashtag. You can find previous tweets that way as shown here. Thank you!


Thank you for your interest in our work!

## Study Design

As a quick overview, the data presented here is from hippocampus (HPC) that spans nine spatial domains plus white matter for a total of sixteen subjects. Each dissection of HPC was designed to center the granular cell layer to well represent the DG. Using this web application you can explore the spatial expression of known genes such as NCDN.

<img src="https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/img/study_design.png" width="1000px" align="left" />

This web application was built such that we could annotate the spots to layers as you can see under the spot-level data tab. Both histologically and gene marker driven manual annotations as well as unsupervised spatial clusters with BayesSpace at k=10 are available. Once we annotated each spot to a layer, we compressed the information by a pseudo-bulking approach into layer-level data minus the layer representing choroid plexus to maximize variance between HPC spatial domains. We then analyzed the expression through a set of models whose results you can also explore through this web application. Finally, you can upload your own gene sets of interest as well as layer enrichment statistics and compare them with our LIBD human lifespan DG Visium dataset.

If you are interested in running this web application locally, you can do so thanks to the spatialLIBD R/Bioconductor package that powers this web application as shown below.

First download the processed spe objects and modeling results here:


## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔭 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)

We provide the following interactive websites, organized by software
labeled with emojis:

- 🔭 `spatialLIBD`
  - [https://libd.shinyapps.io/Lifespan_DG/](https://libd.shinyapps.io/Lifespan_DG/):
    This web application was built such that we could annotate the spots to layers as you can see under the spot-level data tab. Both histologically and gene marker driven manual annotations as well as unsupervised spatial clusters with BayesSpace at k=10 are available. Once we annotated each spot to a layer, we compressed the information by a pseudo-bulking approach into layer-level data minus the layer representing choroid plexus to maximize variance between HPC spatial domains. We then analyzed the expression through a set of models whose results you can also explore through this web application.

### Local `spatialLIBD` apps

If you are interested in running the
[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) applications
locally, you can do so thanks to the
[`spatialLIBD::run_app()`](http://research.libd.org/spatialLIBD/reference/run_app.html),
which you can also use with your own data as shown in our [vignette for
publicly available datasets provided by 10x
Genomics](http://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html).

``` r
## Run this web application locally with:
## Deploy the website
###spatialLIBD::run_app(
###    spe,
###    sce_layer = spe_pseudo,
###    modeling_results = modeling_results,
###    sig_genes = NULL,
###    title = "spatial_DG_lifespan, Visium",
###    spe_discrete_vars = c("BayesSpace", "ManualAnnotation"),
###    spe_continuous_vars = c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio"),
###    default_cluster = "BayesSpace"
###)
## You will have more control about the length of the
## session and memory usage.

## You could also use this function to visualize your
## own data given some requirements described
## in detail in the package vignette documentation
## at http://research.libd.org/spatialLIBD/.

```

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatial_DG_lifespan/issues](https://github.com/LieberInstitute/spatial_DG_lifespan/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Citing our work

Please cite this manuscript
if you use data from this project.


Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.


### Cite `spatialLIBD`

Below is the citation output from using `citation('spatialLIBD')` in R.
Please run this yourself to check for any updates on how to cite
**spatialLIBD**.

``` r
print(citation("spatialLIBD")[1], bibtex = TRUE)
#> 
#> Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K,
#> Maynard KR, Collado-Torres L (2022). "spatialLIBD: an R/Bioconductor
#> package to visualize spatially-resolved transcriptomics data." _BMC
#> Genomics_. doi:10.1186/s12864-022-08601-w
#> <https://doi.org/10.1186/s12864-022-08601-w>,
#> <https://doi.org/10.1186/s12864-022-08601-w>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {spatialLIBD: an R/Bioconductor package to visualize spatially-resolved transcriptomics data},
#>     author = {Brenda Pardo and Abby Spangler and Lukas M. Weber and Stephanie C. Hicks and Andrew E. Jaffe and Keri Martinowich and Kristen R. Maynard and Leonardo Collado-Torres},
#>     year = {2022},
#>     journal = {BMC Genomics},
#>     doi = {10.1186/s12864-022-08601-w},
#>     url = {https://doi.org/10.1186/s12864-022-08601-w},
#>   }
```

Please note that the `spatialLIBD` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing the package.

### Cite `VistoSeg`

To cite [`VistoSeg`](http://research.libd.org/VistoSeg/) please use:

> VistoSeg: processing utilities for high-resolution Visium/Visium-IF
> images for spatial transcriptomics data. Madhavi Tippani, Heena R.
> Divecha, Joseph L. Catallini II, Sang Ho Kwon, Lukas M. Weber, Abby
> Spangler, Andrew E. Jaffe, Stephanie C. Hicks, Keri Martinowich,
> Leonardo Collado-Torres, Stephanie C. Page, Kristen R. Maynard bioRxiv
> 2021.08.04.452489; doi: <https://doi.org/10.1101/2021.08.04.452489>

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article {Tippani2021.08.04.452489,
        author = {Tippani, Madhavi and Divecha, Heena R. and Catallini, Joseph L. and Kwon, Sang Ho and Weber, Lukas M. and Spangler, Abby and Jaffe, Andrew E. and Hicks, Stephanie C. and Martinowich, Keri and Collado-Torres, Leonardo and Page, Stephanie C. and Maynard, Kristen R.},
        title = {VistoSeg: processing utilities for high-resolution Visium/Visium-IF images for spatial transcriptomics data},
        elocation-id = {2021.08.04.452489},
        year = {2022},
        doi = {10.1101/2021.08.04.452489},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2022/05/13/2021.08.04.452489},
        eprint = {https://www.biorxiv.org/content/early/2022/05/13/2021.08.04.452489.full.pdf},
        journal = {bioRxiv}
    }

## Data Access

We highly value open data sharing and believe that doing so accelerates
science, as was the case between our
[`HumanPilot`](https://doi.org/10.1038/s41593-020-00787-0) and the
external [`BayesSpace`](https://doi.org/10.1038/s41587-021-00935-2)
projects, documented [on this
slide](https://speakerdeck.com/lcolladotor/hca-la-2022?slide=18).

### Processed Data

[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) also allows
you to access the data from this project as ready to use R objects. That
is, a:

- [`SpatialExperiment`](https://doi.org/10.1093/bioinformatics/btac299)
  object for the Visium-SPG samples (n = 10)

You can use the
[`zellkonverter`](https://bioconductor.org/packages/zellkonverter/)
Bioconductor package to convert any of them into Python
[`AnnData`](https://anndata.readthedocs.io/en/latest/) objects. If you
browse our code, you can find examples of such conversions.

If you are unfamiliar with these tools, you might want to check the
[LIBD rstats club](http://research.libd.org/rstatsclub/#.Y4hWlOzMJUM)
(check and search keywords on the
[schedule](https://docs.google.com/spreadsheets/d/1is8dZSd0FZ9Qi1Zvq1uRhm-P1McnJRd_zxdAfCRoMfA/edit?usp=sharing))
videos and resources.

#### Installing spatialLIBD

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `spatialLIBD` from
[Bioconductor](http://bioconductor.org/) with the following code:

``` r
## Install BiocManager in order to install Bioconductor packages properly
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
## Check that you have a valid R/Bioconductor installation
BiocManager::valid()
## Now install spatialLIBD from Bioconductor
## (this version has been tested on macOS, winOS, linux)
BiocManager::install("spatialLIBD")
## If you need the development version from GitHub you can use the following:
# BiocManager::install("LieberInstitute/spatialLIBD")
## Note that this version might include changes that have not been tested
## properly on all operating systems.
```

### R objects

Using `spatialLIBD` you can access the spatialDLPFC transcriptomics data
from the 10x Genomics Visium platform. For example, this is the code you
can use to access the spatially-resolved data. For more details, check
the help file for `fetch_data()`.

``` r
## Check that you have a recent version of spatialLIBD installed
stopifnot(packageVersion("spatialLIBD") >= "1.11.12")
## Download the spot-level data
spe <- spatialLIBD::fetch_data(type = "Visium_SPG_AD_Visium_wholegenome_spe")

## This is a SpatialExperiment object
spe
#> class: SpatialExperiment 
#> dim: 27853 38115 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(27853): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
#>   ENSG00000277196
#> rowData names(7): source type ... gene_type gene_search
#> colnames(38115): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1
#>   TTGTTTGTGTAAATTC-1
#> colData names(113): key sample_id ... APOe path_groups_colors
#> reducedDimNames(15): 10x_pca 10x_tsne ... TSNE_perplexity50.HARMONY
#>   TSNE_perplexity80.HARMONY
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
lobstr::obj_size(spe)
#> 2.29 GB

## Remake the logo image
p_pathology <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "path_groups",
    sampleid = "V10A27106_D1_Br3880",
    colors = spe$path_groups_colors[!duplicated(spe$path_groups_colors)],
    spatial = FALSE,
    ... = " Visium SPG AD\nPathology groups -- made with spatialLIBD"
)
p_pathology
```

### Raw data

You can access all the raw data through
[Zenodo](https://doi.org/10.5281/zenodo.10126688)
This includes all the input FASTQ files, as well as the raw images and processed spe objects.

## Internal

- JHPCE locations:
  `/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/`
- Slack channel:
  [`libd_lifespan_hpc_spatial`](https://jhu-genomics.slack.com/archives/C02TV4F1MCM).

### Files:

- `code`: R, python, and shell scripts for running various analyses.
- `plots`: plots generated by R analysis scripts in `.pdf` or `.png`
  format
- `processed-data`
  - `Images`: images used for running `SpaceRanger` and other images
  - `spaceranger`: `SpaceRanger` output files
- `raw-data`
  - `FASTQ`: FASTQ files.
  - `Images`: raw images from the scanner in `.tif` format for each
   slide (around 8GB each). Each slide contains the image
    for four capture areas.

This GitHub repository is organized along the [*R/Bioconductor-powered
Team Data Science* group
guidelines](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.Yaf9fPHMIdk).
It follows the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
structure.

### Other related files

- Reference transcriptome from 10x Genomics:
  `/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/`

