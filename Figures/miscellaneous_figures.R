# Code for other miscellaneous figures

# Dot plot of DMGs from MCI v CU and LOAD vs MCI
library(clusterProfiler)
library(viridis)

load("MCI Figures/MCI CONTROL/Functional Figures/goData.rdata")
dm.genes.df.mci <- dm.genes.df
dm.genes.df.hyper.mci <- dm.genes.df.hyper
dm.genes.df.hypo.mci <- dm.genes.df.hypo
load("MCI Figures/LOAD MCI/Functional_Figures/goData.rdata")
dm.genes.df.hyper.load <- dm.genes.df.hyper
dm.genes.df.hypo.load <- dm.genes.df.hypo
dm.genes.df.load <- dm.genes.df

genes <- list(dm.genes.df.hyper.mci$gene_name, dm.genes.df.hypo.mci$gene_name, dm.genes.df.hyper.load$gene_name, dm.genes.df.hypo.load$gene_name)
names(genes) <- c("MCI - CU (Hyper)", "MCI - CU (Hypo)", "AD - MCI (Hyper)", "AD - MCI (Hypo)")
ck <- compareCluster(genes,fun='enrichGO',OrgDb="org.Hs.eg.db",ont="BP",keyType="SYMBOL")
pdf("go.pdf",height = 15, width = 11)
dotplot(ck, by = "Count", showCategory=10) + scale_color_viridis(option = "turbo")
dev.off()



# Venn diagram
library(data.table)
library(GenomicRanges)

pvals.mci.cu <- fread("~/Desktop/MCI Figures/MCI CONTROL/Functional Figures/pvals.bed", header=T)
pvals.load.cu <- fread("~/Desktop/MCI Figures/LOAD CONTROL/Functional_Figures/pvals.bed", header=T)
pvals.load.mci <- fread("~/Desktop/MCI Figures/LOAD MCI/Functional_Figures/pvals.bed", header=T)

cutoff <- 0.025
pvals.df <- pvals.mci.cu
sig.mci <- pvals.df[which(pvals.df$lfdr.from.ss<0.05 & abs(pvals.df$pi.diff.mci.ctrl)>cutoff),]
pvals.df <- pvals.load.cu
sig.load.cu <- pvals.df[which(pvals.df$lfdr.from.ss<0.05 & abs(pvals.df$pi.diff.load.ctrl)>cutoff),]
pvals.df <- pvals.load.mci
sig.load.mci <- pvals.df[which(pvals.df$lfdr.from.ss<0.05 & abs(pvals.df$pi.diff.load.mci)>cutoff),]

sig.mci.gr <- with(sig.mci, GRanges(chr, IRanges(start, start)))
sig.load.cu.gr <- with(sig.load.cu, GRanges(chr, IRanges(start, start)))
sig.load.mci.gr <- with(sig.load.mci, GRanges(chr, IRanges(start, start)))

overlap.mci.loadCu <- as.data.frame(findOverlaps(sig.mci.gr, sig.load.cu.gr)) 
# 814
overlap.mci.loadMci <- as.data.frame(findOverlaps(sig.mci.gr, sig.load.mci.gr)) 
# 1545
overlap.loadCu.loadMci <- as.data.frame(findOverlaps(sig.load.cu.gr, sig.load.mci.gr)) 
# 1862

sig.mci.gr.subset <- sig.mci.gr[overlap.mci.loadCu[,1],]
overlap.triple <- as.data.frame(findOverlaps(sig.mci.gr.subset, sig.load.mci.gr))
# 58

mci.x <- sig.mci.gr.subset[overlap.triple[,1],]
mci.x <- as.data.frame(mci.x)
runningTotal <- c()
for (i in 1:nrow(mci.x)) {
# MCI vs CU leg
xx <- sig.mci[which(sig.mci$chr == mci.x[i, "seqnames"] & sig.mci$start == mci.x[i,"start"]),]
if(xx$pi.diff.mci.ctrl < 0) {
type1 <- "hypo"}
if(xx$pi.diff.mci.ctrl > 0) {
type1 <- "hyper"}
rm(xx)
# LOAD vs CU leg
xx <- sig.load.cu[which(sig.load.cu$chr == mci.x[i, "seqnames"] & sig.load.cu$start == mci.x[i,"start"]),]
if(xx$pi.diff.load.ctrl < 0) {
type2 <- "hypo"}
if(xx$pi.diff.load.ctrl > 0) {
type2 <- "hyper"}
rm(xx)
# LOAD vs MCI leg
xx <- sig.load.mci[which(sig.load.mci$chr == mci.x[i, "seqnames"] & sig.load.mci$start == mci.x[i,"start"]),]
if(xx$pi.diff.load.mci < 0) {
type3 <- "hypo"}
if(xx$pi.diff.load.mci > 0) {
type3 <- "hyper"}
rm(xx)
xx <- cbind(type1, type2, type3)
runningTotal <- rbind(runningTotal, xx)}
runningTotal <- as.data.frame(runningTotal)
write.table(runningTotal, file = "~/Desktop/totals.txt",quote=F,row.names=F,sep='\t')


sig.overlaps <- sig.mci.gr.subset[overlap.triple[,1],]
sig.overlaps$pi.diff <- 0.025
# generate enhancers.to.test etc in other code
combine_dmps_with_intervals(sig.overlaps, enhancers.to.test)
combine_dmps_with_intervals(sig.overlaps, promoters.to.test)
combine_dmps_with_intervals(sig.overlaps, baits.to.test)

# overlap with transitional DMPs
rm(pvals.mci.cu, pvals.load.mci, pvals.load.cu, pvals.df)
pvals.transitional <- fread("~/Desktop/MCI Figures/CONTROL MCI LOAD/pvals.bed", header=T)
sig.trans <- pvals.transitional[which(pvals.transitional$lfdr.from.ss<0.05 & ((pvals.transitional$pi.diff.mci.ctrl > 0.025 & pvals.transitional$pi.diff.load.mci > 0.025) | (pvals.transitional$pi.diff.mci.ctrl < -0.025 & pvals.transitional$pi.diff.load.mci < -0.025))),]
sig.trans.gr <- with(sig.trans, GRanges(chr, IRanges(start,end)))
x <- as.data.frame(findOverlaps(sig.trans.gr, sig.overlaps))
sig.x <- sig.trans[x[,1],]

# HOMER analysis of 59 overlapping CpGs
#findMotifsGenome.pl overlapping.bed hg38 out -cpg -size given -p 50 -nomotif

# Demographics table tests

# load packages
suppressPackageStartupMessages({
library(data.table)
library(devtools)
library(tidyverse)
library(harmonicmeanp)
library(openxlsx)
library(EnsDb.Hsapiens.v86)
library(ggsci)
library(dplyr)
library(magrittr)
library(webshot)
library(networkD3)
library(GenomicRanges)
})

# read in functions
source_url("https://raw.githubusercontent.com/andymadrid/WGMS_New_Figures_Code/main/functions_all_MCI_v_CTRL.R")

# set up some cutoffs
ALPHA <- 0.05
DMALPHA <- 0.025
DMGENE.ALPHA <- 0.01

master.full.df <- read_csv("WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv", show_col_types = F)
sample.ids <- get_sample_ids_from_dss_inputs(file.path("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/MCI_CONTROL/test-diagnostic-group-coded/","DSS-outputs-chr22.RData"))
master.df <- munge_master_df(master.full.df, sample.ids)                                                         
load.sample.ids <- get_group_ids(master.df, "LOAD")
control.sample.ids <- get_group_ids(master.df, "CONTROL")
mci.sample.ids <- get_group_ids(master.full.df, "MCI")

ancestry.test <- tabulate_and_test(master.df, "race_primary", "fisher")
source.test <- tabulate_and_test(master.df, "source", "chi")
sex.test <- tabulate_and_test(master.df, "sex", "chi")
APOE.test <- tabulate_and_test(master.df, "APOE_risk_allele", "chi")
age.test <- test_continuous_three_group(master.df, "age_at_visit")
bmi.test <- test_continuous_three_group(master.df, "bmi")
education.test <- test_continuous_three_group(master.df, "education")
save(master.df,sample.ids,load.sample.ids,control.sample.ids,mci.sample.ids,ancestry.test,source.test,sex.test,APOE.test,age.test,bmi.test,education.test,file="/media/Data/demographics_test_results.rdata")



# PCA plots for batch effects
library(ggplot2)

load("/media/Data/WGBS/LOAD_MCI/Inputs/CONTROL_MCI_LOAD_New/filtered-DSS-inputs-chr1.RData")
x <- as.data.frame(pca.out$x)
var_explained <- pca.out$sdev^2/sum(pca.out$sdev^2)

pdf("/media/Data/pca.pdf")

# diagnostic group
ggplot(data=x,aes(x=PC1,y=PC2,color=df$diagnostic_group)) +
	geom_point(size=1) +
	theme_classic() +
	theme(text=element_text(size=20)) +
	labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
	scale_color_manual(values=c("firebrick2","darkgoldenrod2","mediumblue")) +
	labs(color = "Group")

# source
ggplot(data=x,aes(x=PC1,y=PC2,color=df$source)) +
	geom_point(size=1) +
	theme_classic() +
	theme(text=element_text(size=20)) +
	labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
	scale_color_manual(values=c("firebrick2","mediumblue")) +
	labs(color = "Source")

# sex
ggplot(data=x,aes(x=PC1,y=PC2,color=df$sex)) +
	geom_point(size=1) +
	theme_classic() +
	theme(text=element_text(size=20)) +
	labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
	scale_color_manual(values=c("firebrick2","mediumblue")) +
	labs(color = "Sex")

# ancestry
ggplot(data=x,aes(x=PC1,y=PC2,color=df$race_primary)) +
	geom_point(size=1) +
	theme_classic() +
	theme(text=element_text(size=20)) +
	labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
	scale_color_manual(values=c("firebrick2","darkgoldenrod2","mediumblue")) +
	labs(color = "Ancestry")

dev.off()
