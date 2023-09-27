########################################################################
# spatial_DG_lifespan project
# DE analysis of pseudo-bulked BayesSpace clusters with age t-statistics
# Anthony Ramnauth, Sept 25 2023
########################################################################

setwd("/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(spatialLIBD)
    library(rafalib)
    library(limma)
    library(SummarizedExperiment)
    library(jaffelab)
    library(sva)
    library(limma)
    library(edgeR)
    library(RColorBrewer)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data", "pseudobulk_spe", "pseudobulk_spe.rds"))

# adapted from Jaffe AE. et al., Nat. Neuro. https://doi.org/10.1038/s41593-020-0604-z
# https://github.com/LieberInstitute/dg_hippo_paper/blob/master/03_run_ageAndDx_analyses_lme.R

#####################
## get qSVs #########
#####################

mod = model.matrix(~ 0 + BayesSpace + age + sex , data = colData(spe_pseudo))

# do qSVA
k = num.sv(assays(spe_pseudo)$logcounts, mod)
qSVs = prcomp(t(assays(spe_pseudo)$logcounts))$x[,1:k]
colnames(qSVs) = paste0("qSV", 1:k)

###############################
######### analysis ############
###############################

###### GENE LEVEL ###########

#### whole HPC #####
mod_HPC = model.matrix(~ age + sex + qSVs,
				data = colData(spe_pseudo))

dge_HPC <- DGEList(counts = assays(spe_pseudo)$counts)
dge_HPC <- calcNormFactors(dge_HPC)
vGene_HPC <- voom(dge_HPC, mod_HPC, plot = TRUE)
geneFit_HPC = lmFit(vGene_HPC, mod_HPC)
geneFit_HPC = eBayes(geneFit_HPC)
geneAgeStats_HPC = topTable(geneFit_HPC, coef = 2, n = nrow(dge_HPC), sort="none")

#############################################################################################

#### GCL only #####
spe_pseudo_GCL <- spe_pseudo[, which(spe_pseudo$BayesSpace == "7")]

#####################
## get qSVs #########
#####################

mod_GCL = model.matrix(~ age + sex , data = colData(spe_pseudo_GCL))

# do qSVA
k_GCL = num.sv(assays(spe_pseudo_GCL)$logcounts, mod_GCL)
qSVs_GCL = prcomp(t(assays(spe_pseudo_GCL)$logcounts))$x[,1:k_GCL]
colnames(qSVs_GCL) = paste0("qSV", 1:k_GCL)

mod_GCL2 = model.matrix(~ age + sex + qSVs_GCL,
				data = colData(spe_pseudo_GCL))

dge_GCL <- DGEList(counts = assays(spe_pseudo_GCL)$counts)
dge_GCL <- calcNormFactors(dge_GCL)
vGene_GCL <- voom(dge_GCL, mod_GCL2, plot = TRUE)
geneFit_GCL = lmFit(vGene_GCL, mod_GCL2)
geneFit_GCL = eBayes(geneFit_GCL)
geneAgeStats_GCL = topTable(geneFit_GCL, coef = 2, n = nrow(dge_GCL), sort="none")

plot(geneAgeStats_GCL$t, geneAgeStats_HPC$t)

#############################################################################################

#### CA&4 only #####
spe_pseudo_CA3_4 <- spe_pseudo[, which(spe_pseudo$BayesSpace == "4")]

#####################
## get qSVs #########
#####################

mod_CA3_4 = model.matrix(~ age + sex , data = colData(spe_pseudo_CA3_4))

# do qSVA
k_CA3_4 = num.sv(assays(spe_pseudo_CA3_4)$logcounts, mod_CA3_4)
qSVs_CA3_4 = prcomp(t(assays(spe_pseudo_CA3_4)$logcounts))$x[,1:k_CA3_4]
colnames(qSVs_CA3_4) = paste0("qSV", 1:k_CA3_4)

mod_CA3_42 = model.matrix(~ age + sex + qSVs_CA3_4,
				data = colData(spe_pseudo_CA3_4))

dge_CA3_4 <- DGEList(counts = assays(spe_pseudo_CA3_4)$counts)
dge_CA3_4 <- calcNormFactors(dge_CA3_4)
vGene_CA3_4 <- voom(dge_CA3_4, mod_CA3_42, plot = TRUE)
geneFit_CA3_4 = lmFit(vGene_CA3_4, mod_CA3_42)
geneFit_CA3_4 = eBayes(geneFit_CA3_4)
geneAgeStats_CA3_4 = topTable(geneFit_CA3_4, coef = 2, n = nrow(dge_CA3_4), sort="none")

plot(geneAgeStats_CA3_4$t, geneAgeStats_HPC$t)

##############################################################################################

#### ML only #####
spe_pseudo_ML <- spe_pseudo[, which(spe_pseudo$BayesSpace == "2")]

#####################
## get qSVs #########
#####################

mod_ML = model.matrix(~ age + sex , data = colData(spe_pseudo_ML))

# do qSVA
k_ML = num.sv(assays(spe_pseudo_ML)$logcounts, mod_ML)
qSVs_ML = prcomp(t(assays(spe_pseudo_ML)$logcounts))$x[,1:k_ML]
colnames(qSVs_ML) = paste0("qSV", 1:k_ML)

mod_ML2 = model.matrix(~ age + sex + qSVs_ML,
				data = colData(spe_pseudo_ML))

dge_ML <- DGEList(counts = assays(spe_pseudo_ML)$counts)
dge_ML <- calcNormFactors(dge_ML)
vGene_ML <- voom(dge_ML, mod_ML2, plot = TRUE)
geneFit_ML = lmFit(vGene_ML, mod_ML2)
geneFit_ML = eBayes(geneFit_ML)
geneAgeStats_ML = topTable(geneFit_ML, coef = 2, n = nrow(dge_ML), sort="none")

plot(geneAgeStats_ML$t, geneAgeStats_HPC$t)

###############################################################################################

#### SGZ only #####
spe_pseudo_SGZ <- spe_pseudo[, which(spe_pseudo$BayesSpace == "6")]

#####################
## get qSVs #########
#####################

mod_SGZ = model.matrix(~ age + sex , data = colData(spe_pseudo_SGZ))

# do qSVA
k_SGZ = num.sv(assays(spe_pseudo_SGZ)$logcounts, mod_SGZ)
qSVs_SGZ = prcomp(t(assays(spe_pseudo_SGZ)$logcounts))$x[,1:k_SGZ]
colnames(qSVs_SGZ) = paste0("qSV", 1:k_SGZ)

mod_SGZ2 = model.matrix(~ age + sex + qSVs_SGZ,
				data = colData(spe_pseudo_SGZ))

dge_SGZ <- DGEList(counts = assays(spe_pseudo_SGZ)$counts)
dge_SGZ <- calcNormFactors(dge_SGZ)
vGene_SGZ <- voom(dge_SGZ, mod_SGZ2, plot = TRUE)
geneFit_SGZ = lmFit(vGene_SGZ, mod_SGZ2)
geneFit_SGZ = eBayes(geneFit_SGZ)
geneAgeStats_SGZ = topTable(geneFit_SGZ, coef = 2, n = nrow(dge_SGZ), sort="none")

plot(geneAgeStats_SGZ$t, geneAgeStats_HPC$t)

#############################################################################################

#### SLM only #####
spe_pseudo_SLM <- spe_pseudo[, which(spe_pseudo$BayesSpace == "1")]

#####################
## get qSVs #########
#####################

mod_SLM = model.matrix(~ age + sex , data = colData(spe_pseudo_SLM))

# do qSVA
k_SLM = num.sv(assays(spe_pseudo_SLM)$logcounts, mod_SLM)
qSVs_SLM = prcomp(t(assays(spe_pseudo_SLM)$logcounts))$x[,1:k_SLM]
colnames(qSVs_SLM) = paste0("qSV", 1:k_SLM)

mod_SLM2 = model.matrix(~ age + sex + qSVs_SLM,
				data = colData(spe_pseudo_SLM))

dge_SLM <- DGEList(counts = assays(spe_pseudo_SLM)$counts)
dge_SLM <- calcNormFactors(dge_SLM)
vGene_SLM <- voom(dge_SLM, mod_SLM2, plot = TRUE)
geneFit_SLM = lmFit(vGene_SLM, mod_SLM2)
geneFit_SLM = eBayes(geneFit_SLM)
geneAgeStats_SLM = topTable(geneFit_SLM, coef = 2, n = nrow(dge_SLM), sort="none")

plot(geneAgeStats_SLM$t, geneAgeStats_HPC$t)

#############################################################################################

#### WM only #####
spe_pseudo_WM <- spe_pseudo[, which(spe_pseudo$BayesSpace == "10")]

#####################
## get qSVs #########
#####################

mod_WM = model.matrix(~ age + sex , data = colData(spe_pseudo_WM))

# do qSVA
k_WM = num.sv(assays(spe_pseudo_WM)$logcounts, mod_WM)
qSVs_WM = prcomp(t(assays(spe_pseudo_WM)$logcounts))$x[,1:k_WM]
colnames(qSVs_WM) = paste0("qSV", 1:k_WM)

mod_WM2 = model.matrix(~ age + sex + qSVs_WM,
				data = colData(spe_pseudo_WM))

dge_WM <- DGEList(counts = assays(spe_pseudo_WM)$counts)
dge_WM <- calcNormFactors(dge_WM)
vGene_WM <- voom(dge_WM, mod_WM2, plot = TRUE)
geneFit_WM = lmFit(vGene_WM, mod_WM2)
geneFit_WM = eBayes(geneFit_WM)
geneAgeStats_WM = topTable(geneFit_WM, coef = 2, n = nrow(dge_WM), sort="none")

plot(geneAgeStats_WM$t, geneAgeStats_HPC$t)

#############################################################################################

#### SR only #####
spe_pseudo_SR <- spe_pseudo[, which(spe_pseudo$BayesSpace == "5")]

#####################
## get qSVs #########
#####################

mod_SR = model.matrix(~ age + sex , data = colData(spe_pseudo_SR))

# do qSVA
k_SR = num.sv(assays(spe_pseudo_SR)$logcounts, mod_SR)
qSVs_SR = prcomp(t(assays(spe_pseudo_SR)$logcounts))$x[,1:k_SR]
colnames(qSVs_SR) = paste0("qSV", 1:k_SR)

mod_SR2 = model.matrix(~ age + sex + qSVs_SR,
				data = colData(spe_pseudo_SR))

dge_SR <- DGEList(counts = assays(spe_pseudo_SR)$counts)
dge_SR <- calcNormFactors(dge_SR)
vGene_SR <- voom(dge_SR, mod_SR2, plot = TRUE)
geneFit_SR = lmFit(vGene_SR, mod_SR2)
geneFit_SR = eBayes(geneFit_SR)
geneAgeStats_SR = topTable(geneFit_SR, coef = 2, n = nrow(dge_SR), sort="none")

plot(geneAgeStats_SR$t, geneAgeStats_HPC$t)

#############################################################################################

#### SL only #####
spe_pseudo_SL <- spe_pseudo[, which(spe_pseudo$BayesSpace == "8")]

#####################
## get qSVs #########
#####################

mod_SL = model.matrix(~ age + sex , data = colData(spe_pseudo_SL))

# do qSVA
k_SL = num.sv(assays(spe_pseudo_SL)$logcounts, mod_SL)
qSVs_SL = prcomp(t(assays(spe_pseudo_SL)$logcounts))$x[,1:k_SL]
colnames(qSVs_SL) = paste0("qSV", 1:k_SL)

mod_SL2 = model.matrix(~ age + sex + qSVs_SL,
				data = colData(spe_pseudo_SL))

dge_SL <- DGEList(counts = assays(spe_pseudo_SL)$counts)
dge_SL <- calcNormFactors(dge_SL)
vGene_SL <- voom(dge_SL, mod_SL2, plot = TRUE)
geneFit_SL = lmFit(vGene_SL, mod_SL2)
geneFit_SL = eBayes(geneFit_SL)
geneAgeStats_SL = topTable(geneFit_SL, coef = 2, n = nrow(dge_SL), sort="none")

plot(geneAgeStats_SL$t, geneAgeStats_HPC$t)

#############################################################################################

#### SL only #####
spe_pseudo_CA1 <- spe_pseudo[, which(spe_pseudo$BayesSpace == "9")]

#####################
## get qSVs #########
#####################

mod_CA1 = model.matrix(~ age + sex , data = colData(spe_pseudo_CA1))

# do qSVA
k_CA1 = num.sv(assays(spe_pseudo_CA1)$logcounts, mod_CA1)
qSVs_CA1 = prcomp(t(assays(spe_pseudo_CA1)$logcounts))$x[,1:k_CA1]
colnames(qSVs_CA1) = paste0("qSV", 1:k_CA1)

mod_CA12 = model.matrix(~ age + sex + qSVs_CA1,
				data = colData(spe_pseudo_CA1))

dge_CA1 <- DGEList(counts = assays(spe_pseudo_CA1)$counts)
dge_CA1 <- calcNormFactors(dge_CA1)
vGene_CA1 <- voom(dge_CA1, mod_CA12, plot = TRUE)
geneFit_CA1 = lmFit(vGene_CA1, mod_CA12)
geneFit_CA1 = eBayes(geneFit_CA1)
geneAgeStats_CA1 = topTable(geneFit_CA1, coef = 2, n = nrow(dge_CA1), sort="none")

plot(geneAgeStats_CA1$t, geneAgeStats_HPC$t)

#############################################################################################

## interaction model of interaction between spatial domain and aging ##

#### joint interaction #####
mod_joint_age = model.matrix(~0 + age*BayesSpace + sex + qSVs, data = colData(spe_pseudo))

dge_joint <- DGEList(counts = assays(spe_pseudo)$counts)
dge_joint <- calcNormFactors(dge_joint)
vGene_joint <- voom(dge_joint, mod_joint_age, plot = TRUE)

## duplicate correlation
corfit_gene <- duplicateCorrelation(vGene_joint$E, mod_joint_age, block=spe_pseudo$sample_id)

## dup corr
geneFit_joint = lmFit(vGene_joint, mod_joint_age, block=spe_pseudo$sample_id,
        correlation = corfit_gene$consensus.correlation)

geneFit_joint = eBayes(geneFit_joint)
geneAgeStats_int = topTable(geneFit_joint, coef = ncol(mod_joint_age),
	n = nrow(dge_joint), sort="none")

######################################## MERGE #############################################

tmp = cbind(geneAgeStats_HPC, geneAgeStats_SLM, geneAgeStats_ML, geneAgeStats_CA3_4,
    geneAgeStats_SR, geneAgeStats_SGZ,geneAgeStats_GCL, geneAgeStats_SL,
    geneAgeStats_CA1, geneAgeStats_WM, geneAgeStats_int)
colnames(tmp) = paste0(colnames(tmp), "_Age_", rep(c(
    "HPC", "SLM", "ML", "CA3&4", "SR", "SGZ", "GCL", "SL", "CA1", "WM", "Inter"
    ), each=6))
geneAgeStats = rowRanges(spe_pseudo)
mcols(geneAgeStats) = cbind(tmp, mcols(geneAgeStats))

# save .Rdata file
save(geneAgeStats, file = here::here("processed-data", "pseudobulk_spe", "geneLevel_age_interaction.rds"))

# load(here::here("processed-data", "pseudobulk_spe", "geneLevel_age_interaction.rds"))

## write out for CSV
out = mcols(geneAgeStats)

# Phew, checks out!
out = out[,c(69,71,72,1:66)]
out = as.data.frame(out)

# directory to save results
dir_outputs <- here("processed-data", "pseudobulk_spe", "age_stats")

fn_out1 <- file.path(dir_outputs, "ageStats_geneLevel_plus_BayesSpaceInteraction")

# Export summary as .csv file
write.csv(out, fn_out1, row.names = FALSE)

#################################
# Check stats and make some plots
#################################

# HPC stats
sum(geneAgeStats$adj.P.Val_Age_HPC < 0.05)
table(sign(geneAgeStats$t_Age_HPC[geneAgeStats$adj.P.Val_Age_HPC < 0.05]))
2^quantile(abs(geneAgeStats$logFC_Age_HPC[geneAgeStats$adj.P.Val_Age_HPC < 0.05])*10, c(.5, 0.25,0.75))
#     50%      25%      75%
# 1.079487 1.050969 1.127122

#########################################################################################################

# GCL stats
sum(geneAgeStats$adj.P.Val_Age_GCL < 0.05)
# 4
table(sign(geneAgeStats$t_Age_GCL[geneAgeStats$adj.P.Val_Age_GCL < 0.05]))
# 1
# 4
2^quantile(abs(geneAgeStats$logFC_Age_GCL[geneAgeStats$adj.P.Val_Age_GCL < 0.05])*10, c(.5, 0.25,0.75))
#      50%      25%      75%
# 1.442056 1.349549 1.719303

## overlap?
table(geneAgeStats$adj.P.Val_Age_GCL < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("GCL", "HPC"))
#        HPC
# GCL     FALSE  TRUE
#   FALSE 12467   347
#   TRUE      1     3

getOR(table(geneAgeStats$adj.P.Val_Age_GCL < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("GCL", "HPC")))
# [1] 107.7839

## same direction?

eitherIndex = which(geneAgeStats$adj.P.Val_Age_GCL < 0.05 | geneAgeStats$adj.P.Val_Age_HPC < 0.05)
bothIndex = which(geneAgeStats$adj.P.Val_Age_GCL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05)
length(eitherIndex)
# [1] 351
tt = table(sign(geneAgeStats$t_Age_GCL[eitherIndex]),sign(geneAgeStats$t_Age_HPC[eitherIndex]) )
#       -1   1
#   -1 126   7
#   1   12 206
sum(diag(tt))
# [1] 332
sum(diag(tt))/sum(tt)
# [1] 0.9458689
(sum(tt)-sum(diag(tt)))
# [1] 19
(sum(tt)-sum(diag(tt)))/sum(tt)
# [1] 0.05413105
ttBoth = table(sign(geneAgeStats$t_Age_GCL[bothIndex]),sign(geneAgeStats$t_Age_HPC[bothIndex]) )
#    1
#  1 3

#########################################################################################################

# WM stats
sum(geneAgeStats$adj.P.Val_Age_WM < 0.05)
# 0
table(sign(geneAgeStats$t_Age_WM[geneAgeStats$adj.P.Val_Age_WM < 0.05]))
2^quantile(abs(geneAgeStats$logFC_Age_WM[geneAgeStats$adj.P.Val_Age_WM < 0.05])*10, c(.5, 0.25,0.75))

## overlap?
table(geneAgeStats$adj.P.Val_Age_WM < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("WM", "HPC"))
#       HPC
# WM      FALSE  TRUE
#   FALSE 12468   350

#########################################################################################################

# SLM stats
sum(geneAgeStats$adj.P.Val_Age_SLM < 0.05)
# 13
table(sign(geneAgeStats$t_Age_SLM[geneAgeStats$adj.P.Val_Age_SLM < 0.05]))
# -1  1
#  3 10
2^quantile(abs(geneAgeStats$logFC_Age_SLM[geneAgeStats$adj.P.Val_Age_SLM < 0.05])*10, c(.5, 0.25,0.75))
#      50%      25%      75%
# 1.242205 1.207216 1.311375

## overlap?
table(geneAgeStats$adj.P.Val_Age_SLM < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("SLM", "HPC"))
#        HPC
# SLM     FALSE  TRUE
#   FALSE 12463   342
#   TRUE      5     8

getOR(table(geneAgeStats$adj.P.Val_Age_SLM < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("SLM", "HPC")))
# [1] 58.30643

## same direction?

eitherIndex = which(geneAgeStats$adj.P.Val_Age_SLM < 0.05 | geneAgeStats$adj.P.Val_Age_HPC < 0.05)
bothIndex = which(geneAgeStats$adj.P.Val_Age_SLM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05)
length(eitherIndex)
# [1] 355
tt = table(sign(geneAgeStats$t_Age_SLM[eitherIndex]),sign(geneAgeStats$t_Age_HPC[eitherIndex]) )
#      -1   1
#  -1 123  48
#  1   17 167
sum(diag(tt))
# [1] 290
sum(diag(tt))/sum(tt)
# [1] 0.8169014
(sum(tt)-sum(diag(tt)))
# [1] 65
(sum(tt)-sum(diag(tt)))/sum(tt)
# [1] 0.1830986
ttBoth = table(sign(geneAgeStats$t_Age_SLM[bothIndex]),sign(geneAgeStats$t_Age_HPC[bothIndex]) )
#     -1 1
#  -1  1 0
#  1   0 7

#########################################################################################################

# SGZ stats
sum(geneAgeStats$adj.P.Val_Age_SGZ < 0.05)
# 32
table(sign(geneAgeStats$t_Age_SGZ[geneAgeStats$adj.P.Val_Age_SGZ < 0.05]))
# -1  1
#  4 28
2^quantile(abs(geneAgeStats$logFC_Age_SGZ[geneAgeStats$adj.P.Val_Age_SGZ < 0.05])*10, c(.5, 0.25,0.75))
#      50%      25%      75%
# 1.208634 1.151754 1.245971

## overlap?
table(geneAgeStats$adj.P.Val_Age_SGZ < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("SGZ", "HPC"))
#       HPC
# SGZ     FALSE  TRUE
#   FALSE 12453   333
#   TRUE     15    17

getOR(table(geneAgeStats$adj.P.Val_Age_SGZ < 0.05,
	geneAgeStats$adj.P.Val_Age_HPC < 0.05,
	dnn = c("SGZ", "HPC")))
# [1] 42.38258

## same direction?

eitherIndex = which(geneAgeStats$adj.P.Val_Age_SGZ < 0.05 | geneAgeStats$adj.P.Val_Age_HPC < 0.05)
bothIndex = which(geneAgeStats$adj.P.Val_Age_SGZ < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05)
length(eitherIndex)
# [1] 365
tt = table(sign(geneAgeStats$t_Age_SLM[eitherIndex]),sign(geneAgeStats$t_Age_HPC[eitherIndex]) )
#      -1   1
#  -1 124  48
#  1   19 174
sum(diag(tt))
# [1] 298
sum(diag(tt))/sum(tt)
# [1] 0.8164384
(sum(tt)-sum(diag(tt)))
# [1] 67
(sum(tt)-sum(diag(tt)))/sum(tt)
# [1] 0.1835616
ttBoth = table(sign(geneAgeStats$t_Age_SLM[bothIndex]),sign(geneAgeStats$t_Age_HPC[bothIndex]) )
#     1
#  1 17

####################################################################################################
# PLOTS
####################################################################################################

pdf(file = here::here("plots", "pseudobulked","age_HPC_vs_tstats.pdf"))

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_SLM > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_SLM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_SLM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "SLM"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_SLM, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "SLM (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_SLM, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SLM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "SLM", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_SLM, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "SLM", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_SLM, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SLM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "SLM", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_ML > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_ML < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_ML < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "ML"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_ML, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "ML (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_ML, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_ML < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "ML", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_ML, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "ML", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_ML, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_ML < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "ML", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_SGZ > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_SGZ < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_SGZ < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "SGZ"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_SGZ, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "SGZ (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_SGZ, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SGZ < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.5, pos = 4)
legend("bottomright", c("HPC", "SGZ", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_SGZ, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "SGZ", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_SGZ, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SGZ < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.5, pos = 4)
legend("bottomright", c("HPC", "SGZ", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$`adj.P.Val_Age_CA3&4` > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$`adj.P.Val_Age_CA3&4` < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$`adj.P.Val_Age_CA3&4` < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "CA3&4"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$`t_Age_CA3&4`, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "CA3&4 (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$`t_Age_CA3&4`, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$`adj.P.Val_Age_CA3&4` < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "CA3&4", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$`logFC_Age_CA3&4`, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "CA3&4", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$`logFC_Age_CA3&4`, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$`adj.P.Val_Age_CA3&4` < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "CA3&4", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_SR > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_SR < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_SR < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "SR"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_SR, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "SR (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_SR, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SR < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "SR", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_SR, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "SR", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$t_Age_SR, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SR < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "SR", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_GCL > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_GCL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_GCL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "GCL"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_GCL, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "GCL (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_GCL, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_GCL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "GCL", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_GCL, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "GCL", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_GCL, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_GCL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "GCL", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_SL > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_SL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_SL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "SL"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_SL, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "SL (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_SL, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "SL", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_SL, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "SL", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_SL, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_SL < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "SL", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_CA1 > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_CA1 < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_CA1 < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "CA1"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_CA1, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "CA1 (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_CA1, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_CA1 < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "CA1", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_CA1, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "CA1", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_CA1, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_CA1 < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "CA1", "Both"),
	pch=15,col=2:4,cex=1.5)

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_WM > 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_WM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_WM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HPC"
names(cols)[cols==3] = "WM"
names(cols)[cols==4] = "Both"

par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_WM, geneAgeStats$t_Age_HPC,
	ylim = c(-10, 10), xlim = c(-10, 10), pch = 21,bg=cols,
	xlab = "WM (age t-stat)", ylab = "HPC (age t-stat)",
	main = "T-statistics")
text(geneAgeStats$t_Age_WM, geneAgeStats$t_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_WM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "WM", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_WM, geneAgeStats$logFC_Age_HPC,
	ylim = c(-0.15, 0.15), xlim = c(-0.15, 0.15), pch = 21,
	bg=cols, xlab = "WM", ylab = "HPC", main = "log2 Fold Changes")
text(geneAgeStats$logFC_Age_WM, geneAgeStats$logFC_Age_HPC,
    labels = ifelse((geneAgeStats$adj.P.Val_Age_WM < 0.05 & geneAgeStats$adj.P.Val_Age_HPC > 0.05),
        geneAgeStats$gene_name, NA), cex = 0.7, pos = 4)
legend("bottomright", c("HPC", "WM", "Both"),
	pch=15,col=2:4,cex=1.5)

dev.off()
