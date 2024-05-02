# Enhancer-Promoter plot

setwd("~/Desktop") 

######################################### 
######################################### 
######################################### 
###### Lollipop/Interaction plots
######################################### 
######################################### 
######################################### 

# load packages
suppressPackageStartupMessages({
library(data.table)
library(DMRichR)
library(annotatr)
library(tidyverse)
library(rtracklayer)
library(trackViewer)
library(dplyr)
})

# get annotations of genes and CpG islands (just do this locally once instead of every time you source this code)
hg38 <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_basicgenes')
hg38.cpgs <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_cpgs')
cds <- hg38[grep("exon",hg38$id)]
cgis <- hg38.cpgs[grep("island",hg38.cpgs$id)]

plotLollipops <- function(gene, transcript, pdfFile, cutoff, interactionFile) {

cat("Reading in and setting up enhancer/promoter data for your gene of interest . . .\n")

### a file of five columns including:
### 1: chr
### 2: promoterStart
### 3: promoterEnd
### 4: enhancerStart
### 5: enhancerEnd
### all gathered from S6 dataset

x <- read.table(interactionFile, header=T)

# get promoter/enhancer data
x$promoterMid <- (x$promoterStart + x$promoterEnd)/2
x$enhancerMid <- (x$enhancerStart + x$enhancerEnd)/2

# set dummy variable for plotting arcs
x$frequency <- 15
anchors <- x
fragment <- x
mid <- x
anchors$point <- "anchor"
fragment$point <- "fragment"
mid$point <- "mid"
anchors$coord <- anchors$promoterMid
fragment$coord <- fragment$enhancerMid
mid$coord <- (mid$promoterMid + mid$enhancerMid)/2
mid$frequency <- 20
xx <- rbind(anchors,fragment,mid)
#xx$group <- rep(c(paste0("coord",1:23)),c(3))
xx$group <- rep(c(paste0("coord",1)),c(3))
#ggplot(xx, aes(x = coord, y = frequency, group=group)) +
#	geom_line(stat="smooth",method = "lm",formula = y ~ poly(x, 2),se = FALSE,lineend="round",size=1,alpha=0.8) +
#	theme_linedraw_noframe +
#	geom_segment(aes(x = promoterStart, y = 15, xend = promoterEnd, yend = 15)) +
#	geom_segment(aes(x = enhancerStart, y = 15, xend = enhancerEnd, yend = 15)) +
#	ylim(0,20)
intData <- xx

# filter just for gene/transcript of interest
cat("Filtering for your gene of interest for DMPs. . .\n")
gene.gr <- cds[which(cds$symbol==gene),] # or whatever gene of interest
gene.gr <- gene.gr[which(gene.gr$tx_id==transcript),] # or whatever transcript of interest
start.gr <- with(fragment[1,],GRanges(chr,IRanges(enhancerStart,enhancerEnd)))

### this needs to be changed to promoterStart/end for forward genes
end.gr <- with(anchors[1,],GRanges(chr,IRanges(promoterStart,promoterEnd)))
gene.gr <- c(start.gr,gene.gr)
gene.df <- as.data.frame(gene.gr)
gene.gr2 <- c(end.gr,gene.gr)
gene.df2 <- as.data.frame(gene.gr2)
gene.df2 <- gene.df2[order(gene.df2$start),]

# get orientation for plotting later on
if (gene.df[1,"start"] > gene.df[nrow(gene.df),"end"]) {
	orientation <- "backwards"
}
if (gene.df[1,"start"] < gene.df[nrow(gene.df),"end"]) {
	orientation <- "forwards"
}

# read in MCI vs CU data and process to granges object
cat("MCI vs CU leg . . .\n")
pvals.df <- pvals.mci.cu
sig <- pvals.df[which(pvals.df$lfdr.from.ss<0.05 & abs(pvals.df$pi.diff.mci.ctrl)>cutoff),]
sig$color <- ifelse(sig$pi.diff.mci.ctrl > 0, "firebrick3", "mediumblue")
sig$shape <- "diamond"
notSig <- pvals.df[which(pvals.df$lfdr.from.ss>0.05),]
notSig$color <- "lightgrey"
notSig$shape <- "circle"
notSig$lfdr.from.ss <- 1
set.seed(714)
notSig <- notSig[sample(nrow(notSig), nrow(notSig)*0.01),]
sig.gr <- with(sig,GRanges(chr,IRanges(end,end)))
notSig.gr <- with(notSig,GRanges(chr,IRanges(end,end)))
cpgs.gr <- c(sig.gr, notSig.gr)
pvals.subset <- rbind(sig,notSig)
cpgs <- as.data.frame(cpgs.gr)
gene.df <- gene.df[order(gene.df$start),]
if (orientation == "backwards") {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] &
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] & 
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
}
if (orientation == "forwards") {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] &
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] & 
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
}
lolly.dmps.gr$color <- lolly.dmps$color
lolly.dmps.gr$score <- -log10(lolly.dmps$lfdr.from.ss)
lolly.dmps.gr$label.parameter.gp <- gpar(col="white")
lolly.dmps.gr$dashline.col <- "white"
lolly.dmps.gr$label.parameter.draw <- FALSE
lolly.dmps.gr$shape <- lolly.dmps$shape
lolly.dmps.gr$cex <- ifelse(lolly.dmps.gr$score>-log10(0.05),0.5,0.05)
names(lolly.dmps.gr) <- lolly.dmps$start+1
lolly.dmps.mci.cu.gr <- lolly.dmps.gr

# read in LOAD vs CU data and process to granges object
cat("LOAD vs CU leg . . .\n")
pvals.df <- pvals.load.cu
sig <- pvals.df[which(pvals.df$lfdr.from.ss<0.05 & abs(pvals.df$pi.diff.load.ctrl)>cutoff),]
sig$color <- ifelse(sig$pi.diff.load.ctrl > 0, "firebrick3", "mediumblue")
sig$shape <- "square"
notSig <- pvals.df[which(pvals.df$lfdr.from.ss>0.05),]
notSig$color <- "lightgrey"
notSig$shape <- "circle"
notSig$lfdr.from.ss <- 1
set.seed(415)
notSig <- notSig[sample(nrow(notSig), nrow(notSig)*0.01),]
sig.gr <- with(sig,GRanges(chr,IRanges(end,end)))
notSig.gr <- with(notSig,GRanges(chr,IRanges(end,end)))
cpgs.gr <- c(sig.gr, notSig.gr)
pvals.subset <- rbind(sig,notSig)
cpgs <- as.data.frame(cpgs.gr)
if (orientation == "backwards") {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] &
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] & 
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
}
if (orientation == "forwards") {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] &
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] & 
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
}
lolly.dmps.gr$color <- lolly.dmps$color
lolly.dmps.gr$score <- -log10(lolly.dmps$lfdr.from.ss)
lolly.dmps.gr$label.parameter.gp <- gpar(col="white")
lolly.dmps.gr$dashline.col <- "white"
lolly.dmps.gr$label.parameter.draw <- FALSE
lolly.dmps.gr$shape <- lolly.dmps$shape
lolly.dmps.gr$cex <- ifelse(lolly.dmps.gr$score>-log10(0.05),0.5,0.05)
names(lolly.dmps.gr) <- lolly.dmps$start+1
lolly.dmps.load.cu.gr <- lolly.dmps.gr

# read in LOAD vs MCI data and process to granges object
cat("LOAD vs MCI leg . . .\n")
pvals.df <- pvals.load.mci
sig <- pvals.df[which(pvals.df$lfdr.from.ss<0.05 & abs(pvals.df$pi.diff.load.mci)>cutoff),]
sig$color <- ifelse(sig$pi.diff.load.mci > 0, "firebrick3", "mediumblue")
sig$shape <- "triangle_point_up"
notSig <- pvals.df[which(pvals.df$lfdr.from.ss>0.05),]
notSig$shape <- "circle"
notSig$color <- "lightgrey"
notSig$lfdr.from.ss <- 1
set.seed(608)
notSig <- notSig[sample(nrow(notSig), nrow(notSig)*0.05),]
sig.gr <- with(sig,GRanges(chr,IRanges(end,end)))
notSig.gr <- with(notSig,GRanges(chr,IRanges(end,end)))
cpgs.gr <- c(sig.gr, notSig.gr)
pvals.subset <- rbind(sig,notSig)
cpgs <- as.data.frame(cpgs.gr)
if (orientation == "backwards") {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] &
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] & 
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
}
if (orientation == "forwards") {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] &
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df2[1,"seqnames"]) & cpgs$end>=gene.df2[1,"start"] & 
		cpgs$end<=gene.df2[nrow(gene.df2),"end"]),]
}
lolly.dmps.gr$color <- lolly.dmps$color
lolly.dmps.gr$score <- -log10(lolly.dmps$lfdr.from.ss)
lolly.dmps.gr$label.parameter.gp <- gpar(col="white")
lolly.dmps.gr$dashline.col <- "white"
lolly.dmps.gr$label.parameter.draw <- FALSE
lolly.dmps.gr$shape <- lolly.dmps$shape
lolly.dmps.gr$cex <- ifelse(lolly.dmps.gr$score>-log10(0.05),0.5,0.05)
names(lolly.dmps.gr) <- lolly.dmps$start+1
lolly.dmps.load.mci.gr <- lolly.dmps.gr

cat("Bringing  it all back home. . .\n")
lolly.dmps.gr <- c(lolly.dmps.mci.cu.gr, lolly.dmps.load.cu.gr, lolly.dmps.load.mci.gr)

# reduce if too many CpGs
if (length(lolly.dmps.gr) > 75) {
	cat("Ruh roh . . . that’s a lot of CpGs. Let’s fix that. . .\n")
	sigCGs <- which(lolly.dmps.gr$score>0)
	notSigCGs <- sample(which(lolly.dmps.gr$score==0),0.2*length(which(lolly.dmps.gr$score==0)))
	subsetCGs <- c(sigCGs,notSigCGs)
	lolly.dmps.gr <- lolly.dmps.gr[subsetCGs]
}

# plot it out
cat("Plotting her out . . .\n")
lolly.dmps.gr2 <- lolly.dmps.gr
names(lolly.dmps.gr2) <- 1:length(lolly.dmps.gr2)
lolly.dmps.gr2 <- as.data.frame(lolly.dmps.gr2)

# change lollipop size to make signficant ones more evident
lolly.dmps.gr2$cex <- ifelse(lolly.dmps.gr2$cex== 0.5, 4, 1)

# change shapes depending on which comparison they came from
lolly.dmps.gr2 <-    lolly.dmps.gr2 %>%
	dplyr::mutate(
        	shape2 = 
		case_when(shape == "diamond" ~ 18,
		shape == "triangle_point_up" ~ 17,
		shape == "circle" ~ 16,
		shape == "square" ~ 15))

# change color of lollipop sticks
lolly.dmps.gr2$stickCol <- ifelse(lolly.dmps.gr2$score == 0, "lightgrey", "black")

# get exons for plotting down below
lines_with_exons <- grep("exon",gene.df2$id)

# plot lollipops
yl <- expression( -log [10] (lFDR))
lollipops <- ggplot(lolly.dmps.gr2, aes(x = start, y = score)) +
#	geom_hline(yintercept=-1, color="black", size =1) +
#	geom_segment(x = gene.df2[1,"start"], y = -1, xend = gene.df2[nrow(gene.df2),"end"], yend = -1, color="darkgoldenrod3") +
	geom_segment(x = gene.df2[min(lines_with_exons),"start"], y = -1, xend = gene.df2[max(lines_with_exons),"end"], yend = -1, color="darkgoldenrod3") +
	geom_segment(x = lolly.dmps.gr2$start, y = -1, xend = lolly.dmps.gr2$start, yend = lolly.dmps.gr2$score, color = lolly.dmps.gr2$stickCol) +
	geom_point(size=lolly.dmps.gr2$cex, shape = lolly.dmps.gr2$shape2, color = lolly.dmps.gr2$color) +
	xlim(gene.df2[1,"start"]-1000, gene.df2[nrow(gene.df2),"end"]+1000) +
	theme_classic() +
	ylab(yl) +
	ylim(-2,20)

# add in exons
for (i in 1:nrow(gene.df)) {
	if (orientation == "forwards") {
		if (i == 1) {
			lollipops <- lollipops +
#			geom_segment(x = gene.df[i,"start"], y = -1, xend = gene.df[i,"end"], yend = -1, color = "black", size = 0.5)
#			geom_rect(xmin = gene.df[i,"start"], ymin = -1, xmax = gene.df[i,"end"], ymax = -1, color = "black", fill = "black")
			cat("Meow\n")		
		}
		else {
			lollipops <- lollipops +
#			geom_segment(x = gene.df[i,"start"], y = -1.5, xend = gene.df[i,"end"], yend = -0.5, color = "darkgoldenrod2", size = 0.5)
			geom_rect(xmin = gene.df[i,"start"], ymin = -1.5, xmax = gene.df[i,"end"], ymax = -0.5, color = "darkgoldenrod2", fill = "darkgoldenrod3")

		}
	}
	if (orientation == "backwards") {
		if (i != nrow(gene.df)) {
			lollipops <- lollipops +
#			geom_segment(x = gene.df[i,"start"], y = -1.5, xend = gene.df[i,"end"], yend = -0.5, color = "darkgoldenrod2", size = 0.5)
			geom_rect(xmin = gene.df[i,"start"], ymin = -1.5, xmax = gene.df[i,"end"], ymax = -0.5, color = "darkgoldenrod2", fill = "darkgoldenrod2")

		}
		else {
			lollipops <- lollipops +
#			geom_segment(x = gene.df[i,"start"], y = -1, xend = gene.df[i,"end"], yend = -1, color = "black", size = 0.5)
#			geom_rect(xmin = gene.df[i,"start"], ymin = -1, xmax = gene.df[i,"end"], ymax = -1, color = "black", fill = "black")
			cat("Meow\n")		

		}
	}
}
# add in enhancer-promoter interactions
lollipops <- lollipops + geom_segment(x = xx[1,"promoterStart"], y = 10, xend = xx[1,"promoterEnd"], yend = 10)
lollipops <- lollipops + geom_segment(x = xx[1,"enhancerStart"], y = 10, xend = xx[1,"enhancerEnd"], yend = 10)

# add in enhancer-promoter interaction curve
if (orientation == "forwards") {
	lollipops <- lollipops + geom_curve(x = xx[1,"coord"], y = 10, xend = xx[2,"coord"], yend = 10, ncp=3000)
}
if (orientation == "backwards") {
	lollipops <- lollipops + geom_curve(x = xx[2,"coord"], y = 10, xend = xx[1,"coord"], yend = 10, ncp=3000)
}

# add chromosome to plot
lollipops <- lollipops + xlab(xx[1,"chr"])


# add CpG islands (if they’re there)
cat("Time to get tropical . . .\n")
cgis.overlap <- as.data.frame(findOverlaps(cgis,gene.gr))
if (nrow(cgis.overlap) > 0) {
	cgis.sub <- as.data.frame(cgis[cgis.overlap[,1],])
	for (i in 1:nrow(cgis.sub)) {
		lollipops <- lollipops +
			geom_rect(xmin = cgis.sub[i,"start"], ymin = -2, xmax = cgis.sub[i,"end"], ymax = -1.6, color = "forestgreen", fill = "forestgreen")
	}
}
if (nrow(cgis.overlap) == 0) {
	cat("\tJay kay, bro, you ain’t got any CGIs here . . .\n")
}


# find point to add in arrows for directionality
gene.df.exons <- gene.df2[lines_with_exons,]
gene.df.exons <- gene.df.exons[order(gene.df.exons$start),]
	
intronArrows <- c()
for (i in 1:(nrow(gene.df.exons)-1)) {
	intronStart <- gene.df.exons[i,"end"]
	intronEnd <- gene.df.exons[i+1,"start"]
	midPoint <- (intronStart+intronEnd)/2
	intronArrows <- rbind(intronArrows, midPoint)
}

# add in arrows to show directionality
for (i in 1:nrow(intronArrows)) {
	if (gene.df.exons[1,"strand"] == "+") {
		lollipops <- lollipops + 
			geom_segment(x = intronArrows[i,1]-1, y = -1, xend = intronArrows[i,1], yend = -1,
			arrow = arrow(length = unit(0.1, "cm"),ends="last"))
	}
	if (gene.df.exons[1,"strand"] == "-") {
		lollipops <- lollipops + 
			geom_segment(x = intronArrows[i,1]-1, y = -1, xend = intronArrows[i,1], yend = -1,
			arrow = arrow(length = unit(0.1, "cm"),ends="first"),color="black")
	}
}

# make things bigger
lollipops <- lollipops + theme(text=element_text(size=20,color="black"))


pdf(pdfFile, height = 8, width = 12)
print(lollipops)
dev.off()

}

pvals.mci.cu <-fread("~/Desktop/MCI Figures/MCI CONTROL/Functional Figures/pvals.bed", header=T)
pvals.load.cu <- fread("~/Desktop/MCI Figures/LOAD CONTROL/Functional_Figures/pvals.bed",header=T)
pvals.load.mci <- fread("~/Desktop/MCI Figures/LOAD MCI/Functional_Figures/pvals.bed",header=T)

# Lollipop plots of AD-associated genes
#plotLollipops("SPINK4", "ENST00000379721.4", "~/Desktop/enhancerInt.SPINK4.pdf", 0.025, "~/Desktop/spink4_coords.txt")
#plotLollipops("B4GALT1", "ENST00000379731.5", "~/Desktop/dmps-lollyplot-B4GALT1.pdf", 0.025)
#plotLollipops("SND1", "ENST00000354725.8", "~/Desktop/enhancerInt.SND1.pdf", 0.025, "~/Desktop/snd1_coords.txt")
#plotLollipops("SNORA70", "ENST00000354725.8", "~/Desktop/enhancerInt.SNORA70.pdf", 0.025, "~/Desktop/snora70_coords.txt")
plotLollipops("TXNIP", "ENST00000582401.6", "~/Desktop/enhancerInt.TXNIP.pdf", 0.025, "~/Desktop/txnip_coords.txt")

