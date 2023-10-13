####################################################################
# Checking Jaffe et al, 2020 Nat. Neuro. for similar gene enrichment
####################################################################

## load libraries
library(jaffelab)
library(SummarizedExperiment)
library(edgeR)
library(recount)
library(genefilter)
library(RColorBrewer)
library(lmerTest)

## make tables
dir.create("plots")
dir.create("tables")

## load data (can be downloaded here: https://research.libd.org/dg_hippo_paper/data.html )
load("merged_dg_hippo_allSamples_n596.rda")

## make factor
colData(rse_gene_joint)$Dataset = ifelse(colData(rse_gene_joint)$Dataset == "DG", "DG-GCL", "HIPPO")
colData(rse_gene_joint)$Dataset = factor(colData(rse_gene_joint)$Dataset)
colData(rse_gene_joint)$Dataset = relevel(colData(rse_gene_joint)$Dataset , "HIPPO")

## filter on RPKM level
geneIndex = rowMeans(getRPKM(rse_gene_joint[,rse_gene_joint$Dataset == "DG-GCL"], "Length")) > 0.5
rse_gene_joint = rse_gene_joint[geneIndex,]

##### DGE ######
dge = DGEList(counts = assays(rse_gene_joint)$counts,
	genes = rowData(rse_gene_joint))
dge = calcNormFactors(dge)

####################################
## start w/ gene level

## modeling
mod = model.matrix(~Dataset + mitoRate + rRNA_rate + overallMapRate + totalAssignedGene ,
	data=colData(rse_gene_joint))

## mean-variance
vGene = voom(dge,mod,plot=FALSE)
gene_dupCorr = duplicateCorrelation(vGene$E, mod, block=rse_gene_joint$BrNum)
save(gene_dupCorr, file = "rdas/geneLevel_duplicateCorrelation_regionVsCell.rda")
# load("rdas/geneLevel_duplicateCorrelation_regionVsCell.rda")

## do analysis
fitGene = lmFit(vGene, block=rse_gene_joint$BrNum,
	correlation=gene_dupCorr$consensus.correlation)

## top table
eBGene = eBayes(fitGene)
outGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene_joint),
	adjust.method = "bonferroni")
outGene$sigColor = as.numeric(outGene$adj.P.Val < 0.01)+1
sigGene = outGene[outGene$adj.P.Val < 0.01,]

# for abstract
dim(sigGene)
max(sigGene$P.Value)

## write out
save(outGene, file = "rdas/DE_output_cellType_lmer.rda")
# load("rdas/DE_output_cellType_lmer.rda")

## subset
theGenes = c("PROX1", "CALB1", "GAD1","GAD2", "SLC17A7","LAMP5", "POSTN", "CD74",
    "DCX", "CCK")
gStats = outGene[match(theGenes, outGene$Symbol),]

## filter
sigGeneDf = as.data.frame(sigGene)
sigGeneDf$gencodeTx = sapply(sigGeneDf$gencodeTx, paste, collapse=";")

write.csv(sigGeneDf[,c(2,5, 9,11:14,10, 1, 3:4,6)], quote=FALSE,
	row.names=FALSE, file="tables/DE_output_cellType_lmer.csv")

####################
##### plots ########
bIndexes = splitit(rse_gene_joint$BrNum)

## pca
pca = prcomp(t(vGene$E))
pcaVars = getPcaVars(pca)

pdf("plots/pcaPlot_enrichment_n224_paired.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(pca$x, pch = 21, bg = rse_gene_joint$Dataset,cex=1.5,
	xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"))
for(j in seq(along=bIndexes)) {
	ii = bIndexes[[j]]
	lines(pca$x[ii,1], pca$x[ii,2], col ="grey",lwd=0.4)
}
# legend("bottomright", levels(rse_gene_joint$Dataset),
	# col=1:2,pch=15,cex=2,nc=2)
dev.off()

### target genes #######
pdf("LAMP5markerGene_enrichment.pdf")
exprs = vGene$E[match(theGenes, vGene$genes$Symbol),]
palette(brewer.pal(5,"Set1"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=3)
for(i in seq(along=theGenes)) {
	boxplot(exprs[i,] ~ rse_gene_joint$Dataset,
		xlab = "", outline=FALSE,	ylab = "Normalized Expression",
		main = theGenes[i], ylim = c(0,10))
	xx = jitter(as.numeric(rse_gene_joint$Dataset),amount=0.1)
	for(j in seq(along=bIndexes)) {
		lines(exprs[i,] ~ xx, data=colData(rse_gene_joint),
			subset=bIndexes[[j]], col ="grey",lwd=0.4)
	}
	points(exprs[i,] ~ xx,	pch = 21, bg = rse_gene_joint$Dataset)
	pv = paste0("p=",signif( gStats$P.Value[i],3))
	legend("bottom", pv, cex=1.6)
}
dev.off()
