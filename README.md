# spatial_DG_lifespan
Spatial transcriptomics (Visium) in dentate gyrus of hippocampus over the lifespan

spatial_DG_lifespan
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/377886452.svg)](https://zenodo.org/badge/latestdoi/377886452)

## Overview

<img src="https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/img/study_overview.png?raw=true" align="left" width="300px" />

Welcome to the Spatial_DG_lifespan project! It is composed of: <br/>
&nbsp;1. a shiny web application that we are hosting at [https://libd.shinyapps.io/Lifespan_DG/](https://libd.shinyapps.io/Lifespan_DG/) that can handle a limited set of concurrent users,  
&nbsp;2. and a research article with the scientific knowledge we drew from this dataset. The analysis code for our project is available [here](https://github.com/LieberInstitute/spatial_DG_lifespan). 

The web application allows you to browse the LIBD human lifespan dentate gyrus (DG) spatial transcriptomics data generated with the 10x Genomics Visium platform. Please check the manuscript or bioRxiv pre-print for more details about this research.
If you tweet about this website, the data or the R package please use the #spatialLIBD hashtag. You can find previous tweets that way as shown here. Thank you!


Thank you for your interest in our work!

<br/>

## Study Design

As a quick overview, the data presented here is from hippocampus (HPC) that spans nine spatial domains plus white matter for a total of sixteen subjects. Each dissection of HPC was designed to center the granular cell layer to well represent the DG. Using this web application you can explore the spatial expression of known genes such as NCDN.

<img src="https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/img/study_design.png?raw=true" width="1000px" align="left" />

This web application was built such that we could annotate the spots to layers as you can see under the spot-level data tab. Both histologically and gene marker driven manual annotations as well as unsupervised spatial clusters with BayesSpace at k=10 are available. Once we annotated each spot to a layer, we compressed the information by a pseudo-bulking approach into layer-level data minus the layer representing choroid plexus to maximize variance between HPC spatial domains. We then analyzed the expression through a set of models whose results you can also explore through this web application. Finally, you can upload your own gene sets of interest as well as layer enrichment statistics and compare them with our LIBD human lifespan DG Visium dataset.

If you are interested in running this web application locally, you can do so thanks to the spatialLIBD R/Bioconductor package that powers this web application as shown below.

First download the processed spe objects and modeling results here:


## Interactive Websites

We provide the following interactive websites, organized by software
labeled with emojis:

- 🔭 `spatial_DG_lifespan`
  - [https://libd.shinyapps.io/Lifespan_DG/](https://libd.shinyapps.io/Lifespan_DG/):
    This web application was built such that we could annotate the spots to layers as you can see under the spot-level data tab. Both histologically and gene marker driven manual annotations as well as unsupervised spatial clusters with BayesSpace at k=10 are available. Once we annotated each spot to a layer, we compressed the information by a pseudo-bulking approach into layer-level data minus the layer representing choroid plexus to maximize variance between HPC spatial domains. We then analyzed the expression through a set of models whose results you can also explore through this web application.

### Local `spatialLIBD` apps

If you are interested in running the
[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) applications
locally, you can do so thanks to the
[`spatialLIBD::run_app()`](http://research.libd.org/spatialLIBD/reference/run_app.html).
First make sure to download the data from [Zenodo](https://doi.org/10.5281/zenodo.10126688).
Then you can uncompress the files and store the spe object in your chosen directory. Then use readRDS(here::here()) to load the spe object in R.
For example:

``` r
## Run this web application locally with:
## spe <- readRDS(here::here("processed-data", "spe.rds"))
## Deploy the website
###spatialLIBD::run_app(
###    spe,
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

### Cell Reports

Anthony D. Ramnauth, Madhavi Tippani, Heena R. Divecha, Alexis R. Papariello, Ryan A. Miller, Elizabeth A. Pattie, Joel E. Kleinman, Kristen R. Maynard, Leonardo Collado-Torres, Thomas M. Hyde, Keri Martinowich, Stephanie C. Hicks, Stephanie C. Page, Spatially-resolved transcriptomics of human dentate gyrus across postnatal lifespan reveals heterogeneity in markers for proliferation, extracellular matrix, and neuroinflammation. Cell Reports 10.1016/j.celrep.2025.115300; doi: http://dx.doi.org/10.1016/j.celrep.2025.115300.

Here's the citation information on [BibTeX](https://en.wikipedia.org/wiki/BibTeX) format.

```
@article{Ramnauth2025,
  title = {Spatiotemporal analysis of gene expression in the human dentate gyrus reveals age-associated changes in cellular maturation and neuroinflammation},
  volume = {44},
  ISSN = {2211-1247},
  url = {http://dx.doi.org/10.1016/j.celrep.2025.115300},
  DOI = {10.1016/j.celrep.2025.115300},
  number = {2},
  journal = {Cell Reports},
  publisher = {Elsevier BV},
  author = {Ramnauth,  Anthony D. and Tippani,  Madhavi and Divecha,  Heena R. and Papariello,  Alexis R. and Miller,  Ryan A. and Nelson,  Erik D. and Thompson,  Jacqueline R. and Pattie,  Elizabeth A. and Kleinman,  Joel E. and Maynard,  Kristen R. and Collado-Torres,  Leonardo and Hyde,  Thomas M. and Martinowich,  Keri and Hicks,  Stephanie C. and Page,  Stephanie C.},
  year = {2025},
  month = feb
}

```


## Data Access

### BioProject
FASTQ data is available at the [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1043093) page for this project. 

### Raw & Processed Data

You can access all raw data, that is the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files, and processed data through [Zenodo](https://doi.org/10.5281/zenodo.10126688).
Processed data includes two spe objects and raw images for the haematoxylin and eosin-stained tissue:
The basic spe object stores the raw counts for all the Visium experiments and metadata.
The processed spe object is after quality control, batch correction, data integration, unsupervised spatial clustering, cell-type deconvolution, and also contains the metadata.

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

