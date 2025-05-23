###########################################################
# spatial_DG_lifespan project
# spot-level deconvolution using cell2location
# Anthony Ramnauth, April 19 2023
###########################################################

# Adapting from Lukas' script https://github.com/lmweber/locus-c/blob/main/code/analyses/07_deconvolution/LC_deconvolution_cell2location_merged.py

# references:
# - cell2location tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
# - Scanpy AnnData conversion from R: https://theislab.github.io/scanpy-in-R/
# - reticulate: https://rstudio.github.io/reticulate/

# cell2location installation (note: from specific git commit):
# https://github.com/BayraktarLab/cell2location


# The following code is intended to be run interactively in a Python REPL
# session running within an R session (see links above). This allows us to
# access previously saved R objects and convert these to Python AnnData objects.

# Lukas used -
# start interactive session on JHPCE GPU node:
# qrsh -l caracol,mem_free=128G,h_vmem=128G

# I'm using -
# qrsh -l bluejay,mem_free=230G,h_vmem=230G
# I might have to use a gpu node as cpu looks like it 2 orders of magnitude slower

# module load conda_R/4.2
# R


# ----------------------
# load Visium SPE object
# ----------------------

# load Visium SPE object for use within Python session with reticulate

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
})

# Load SPE
spe <- readRDS(here::here("processed-data", "harmony_processed_spe", "harmony_spe.rds"))

dim(spe)

# load snRNA-seq SCE object for use within Python session with reticulate

# Load SCE
sce <- readRDS(file = here::here("processed-data", "sce", "sce_sestan_DG_final.rds"))

dim(sce)

# check annotated clusters
table(colData(sce)$Cell_Type)

dim(sce)

# ---------------------------------
# extract components from R objects
# ---------------------------------

# extract components from SPE/SCE objects to create new AnnData objects
# see https://theislab.github.io/scanpy-in-R/, section 4.4.1: "Creating AnnData from SingleCellExperiment"
# note: use standard data.frames instead of DataFrames so reticulate can access them

# SPE object
spe_counts <- assay(spe, "counts")
spe_logcounts <- assay(spe, "logcounts")
spe_rowdata <- as.data.frame(rowData(spe))
spe_coldata <- as.data.frame(colData(spe))
spe_coldata_combined <- cbind(as.data.frame(colData(spe)), spatialCoords(spe))

# SCE object
sce_counts <- assay(sce, "counts")
sce_logcounts <- assay(sce, "logcounts")
sce_rowdata <- as.data.frame(rowData(sce))
sce_coldata <- as.data.frame(colData(sce))

# --------------------
# start Python session
# --------------------

# start interactive Python session with reticulate
# using JHPCE cell2location virtual environment

library(reticulate)
use_condaenv("/jhpce/shared/jhpce/libd/cell2location/0.8a0/cell2location_env", required = TRUE)
reticulate::repl_python()

# ----------------------------
# start cell2location workflow
# ----------------------------

# using Python code adapted from cell2location tutorial
# note requires GPU for faster runtime

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

# define results folders

results_folder = '/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/processed-data/Cell_Type_Deconvolution'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# ----------------------
# create AnnData objects
# ----------------------

# create AnnData objects from R objects with reticulate
# see https://theislab.github.io/scanpy-in-R/, section 4.4.1: "Creating AnnData from SingleCellExperiment"

# note: using raw counts for cell2location

# SPE object
adata_spe = sc.AnnData(
    X = r.spe_counts.T,
    obs = r.spe_coldata_combined,
    var = r.spe_rowdata,
    dtype = r.spe_counts.dtype
)

# SCE object
adata_sce = sc.AnnData(
    X = r.sce_counts.T,
    obs = r.sce_coldata,
    var = r.sce_rowdata,
    dtype = r.sce_counts.dtype
)

# -------------------------------
# continue cell2location workflow
# -------------------------------

# continue cell2location workflow using objects above


# -------------
# preprocessing
# -------------

# update Visium AnnData object to match structure in cell2location tutorial

# rename object
adata_vis = adata_spe
# rename 'sample_id' column
adata_vis.obs['sample'] = adata_vis.obs['sample_id']
# rename 'gene_id' column
adata_vis.var['SYMBOL'] = adata_vis.var['gene_id']

adata_vis.var_names = adata_vis.var['gene_id']
adata_vis.var_names.name = None

# note: mitochondrial genes have already been filtered out; otherwise see
# cell2location tutorial for how to filter them out


# update snRNA-seq AnnData object to match structure in cell2location tutorial

# rename object
adata_ref = adata_sce

# recommended gene filtering from cell2location tutorial

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter object
adata_ref = adata_ref[:, selected].copy()
adata_ref

# --------------------------------------------
# estimation of reference cell type signatures
# --------------------------------------------

# prepare AnnData object for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10x reaction / sample / batch
                        batch_key='sample_ID',
                        # cell type, covariate used for constructing signatures
                        labels_key='Cell_Type'
                       )

# create and train the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

mod.view_anndata_setup()

# use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
#mod.plot_history(20)

# in this section, we export the estimated cell abundance (summary of the posterior distribution)
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# save model
mod.save(f"{ref_run_name}", overwrite=True)


# save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


# examine QC plots
#mod.plot_QC()

# the model and output h5ad can be loaded later like this:
# mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_ref = sc.read_h5ad(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# ------------------------------
# Cell2location: spatial mapping
# ------------------------------

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare AnnData for Cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")


# parameter values:
# N_cells_per_location = 3 for human brain data (maybe more for GCL), instead of default = 30 (tissue cell density)
# detection_alpha = 20 for human brain data

# for more details see tutorial docs:
# https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html


# create and train the model

mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection
    detection_alpha=20
)

mod.view_anndata_setup()

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training'])

# in this section, we export the estimated cell abundance (summary of the posterior distribution)
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

# the model and output h5ad can be loaded later like this:
# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# adata_file = f"{run_name}/sp.h5ad"
# adata_vis = sc.read_h5ad(adata_file)

# examine reconstruction accuracy to assess if there are any issues with mapping
# the plot should be roughly diagonal, strong deviations will signal problems
mod.plot_QC()
mod.plot_spatial_QC_across_batches()
plt.savefig('plot1.png')

# --------------------------
# export AnnData object to R
# --------------------------

# export components of AnnData object back to R to add to SPE object

obs = adata_vis.obs
var = adata_vis.var
uns = adata_vis.uns
obsm = adata_vis.obsm

# close Python REPL session and go back to R session
exit

# now can access Python objects using reticulate 'py$' syntax
str(py$obs)
str(py$var)
names(py$uns)
str(py$obsm)

# main results are stored in 'py$obsm' as follows
str(py$obsm['means_cell_abundance_w_sf'])
str(py$obsm['stds_cell_abundance_w_sf'])
str(py$obsm['q05_cell_abundance_w_sf'])
str(py$obsm['q95_cell_abundance_w_sf'])

# using posterior means
# note: cell2location tutorial uses 5% quantile ('q05') of posterior distribution

# check results
head(py$obsm['means_cell_abundance_w_sf'])
dim(py$obsm['means_cell_abundance_w_sf'])
colnames(py$obsm['means_cell_abundance_w_sf'])
dim(spe)
all(rownames(py$obsm['means_cell_abundance_w_sf']) == colnames(spe))
length(rownames(py$obsm['means_cell_abundance_w_sf']) == colnames(spe))

# cell type abundances (number of cells per spot)
summary(py$obsm['means_cell_abundance_w_sf'])
# deciles
sapply(py$obsm['means_cell_abundance_w_sf'], quantile, seq(0, 1, by = 0.1))

# add cell2location results to SPE object
colData(spe) <- cbind(colData(spe), py$obsm['means_cell_abundance_w_sf'])

# save SPE object for further plotting
saveRDS(spe, file = here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))
