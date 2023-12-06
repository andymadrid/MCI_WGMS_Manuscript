# Make lollipop plots for genes of interest

# load packages
library(data.table)
library(DMRichR)
library(annotatr)

# get annotations
#hg38 <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_basicgenes')
#hg38.cpgs <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_cpgs')
#cds <- hg38[grep("exon",hg38$id)]
#cgis <- hg38.cpgs[grep("island",hg38.cpgs$id)]

plotLollipops <- function(gene, transcript, pdfFile, cutoff) {

# load packages
suppressPackageStartupMessages({
library(trackViewer)
library(dplyr)
library(ggplot2)
})

# filter just for gene/transcript of interest
cat("Filtering for your gene of interest . . .\n")
gene.gr <- cds[which(cds$symbol==gene),] # or whatever gene of interest
gene.gr <- gene.gr[which(gene.gr$tx_id==transcript),] # or whatever transcript of interest
gene.df <- as.data.frame(gene.gr)

# read in MCI vs CU data and process to granges object
cat("MCI leg . . .\n")
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
if (gene.df[1,"start"] < gene.df[nrow(gene.df),"end"]) {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end>=gene.df[1,"start"] &
		cpgs$end<=gene.df[nrow(gene.df),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end>=gene.df[1,"start"] & 
		cpgs$end<=gene.df[nrow(gene.df),"end"]),]
}
if (gene.df[1,"start"] > gene.df[nrow(gene.df),"end"]) {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end<=gene.df[1,"end"] &
		cpgs$end>=gene.df[nrow(gene.df),"start"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end<=gene.df[1,"end"] & 
		cpgs$end>=gene.df[nrow(gene.df),"start"]),]
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
cat("LOAD leg . . .\n")
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
if (gene.df[1,"start"] < gene.df[nrow(gene.df),"end"]) {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end>=gene.df[1,"start"] &
		cpgs$end<=gene.df[nrow(gene.df),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end>=gene.df[1,"start"] & 
		cpgs$end<=gene.df[nrow(gene.df),"end"]),]
}
if (gene.df[1,"start"] > gene.df[nrow(gene.df),"end"]) {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end<=gene.df[1,"end"] &
		cpgs$end>=gene.df[nrow(gene.df),"start"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end<=gene.df[1,"end"] & 
		cpgs$end>=gene.df[nrow(gene.df),"start"]),]
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
if (gene.df[1,"start"] < gene.df[nrow(gene.df),"end"]) {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end>=gene.df[1,"start"] &
		cpgs$end<=gene.df[nrow(gene.df),"end"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end>=gene.df[1,"start"] & 
		cpgs$end<=gene.df[nrow(gene.df),"end"]),]
}
if (gene.df[1,"start"] > gene.df[nrow(gene.df),"end"]) {
	lolly.dmps.gr <- cpgs.gr[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end<=gene.df[1,"end"] &
		cpgs$end>=gene.df[nrow(gene.df),"start"]),]
	lolly.dmps <- pvals.subset[which(cpgs$seqnames==as.character(gene.df[1,"seqnames"]) & cpgs$end<=gene.df[1,"end"] & 
		cpgs$end>=gene.df[nrow(gene.df),"start"]),]
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

# plot lollipops
yl <- expression( -log [10] (lFDR))
if (gene.df[1,"strand"] == "+") {
lollipops <- ggplot(lolly.dmps.gr2, aes(x = start, y = score)) +
#	geom_hline(yintercept=-1, color="darkgoldenrod3", size =1) +
	geom_segment(x = gene.df[1,"start"], y = -1, xend = gene.df[nrow(gene.df),"end"], yend = -1, color="darkgoldenrod3") +
	geom_segment(x = lolly.dmps.gr2$start, y = -1, xend = lolly.dmps.gr2$start, yend = lolly.dmps.gr2$score, color = lolly.dmps.gr2$stickCol) +
	geom_point(size=lolly.dmps.gr2$cex, shape = lolly.dmps.gr2$shape2, color = lolly.dmps.gr2$color) +
	xlim(gene.df[1,"start"]-1000, gene.df[nrow(gene.df),"end"]+1000) +
	theme_classic() +
	ylab(yl) +
	ylim(-2,max(lolly.dmps.gr2$score)+2)
}
if (gene.df[1,"strand"] == "-") {
lollipops <- ggplot(lolly.dmps.gr2, aes(x = start, y = score)) +
#	geom_hline(yintercept=-1, color="darkgoldenrod3", size =1) +
	geom_segment(x = gene.df[nrow(gene.df),"start"], y = -1, xend = gene.df[1,"end"], yend = -1, color="darkgoldenrod3") +
	geom_segment(x = lolly.dmps.gr2$start, y = -1, xend = lolly.dmps.gr2$start, yend = lolly.dmps.gr2$score, color = lolly.dmps.gr2$stickCol) +
	geom_point(size=lolly.dmps.gr2$cex, shape = lolly.dmps.gr2$shape2, color = lolly.dmps.gr2$color) +
	xlim(gene.df[nrow(gene.df),"start"]-1000, gene.df[1,"end"]+1000) +
	theme_classic() +
	ylab(yl) +
	ylim(-2,max(lolly.dmps.gr2$score)+2)
}

# add in exons
for (i in 1:nrow(gene.df)) {
	lollipops <- lollipops +
	geom_rect(xmin = gene.df[i,"start"], ymin = -1.5, xmax = gene.df[i,"end"], ymax = -0.5, color = "darkgoldenrod3", fill = "darkgoldenrod3")
}


# add chromosome to plot
lollipops <- lollipops + xlab(gene.df[1,"seqnames"])

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
intronArrows <- c()
for (i in 1:(nrow(gene.df)-1)) {
	intronStart <- gene.df[i,"end"]
	intronEnd <- gene.df[i+1,"start"]
	midPoint <- (intronStart+intronEnd)/2
	intronArrows <- rbind(intronArrows, midPoint)
}

# add in arrows to show directionality
for (i in 1:nrow(intronArrows)) {
	if (gene.df[1,"strand"] == "+") {
		lollipops <- lollipops + 
			geom_segment(x = intronArrows[i,1]-1, y = -1, xend = intronArrows[i,1], yend = -1,
			arrow = arrow(length = unit(0.5, "cm"),ends="last"))
	}
	if (gene.df[1,"strand"] == "-") {
		lollipops <- lollipops + 
			geom_segment(x = intronArrows[i,1]-1, y = -1, xend = intronArrows[i,1], yend = -1,
			arrow = arrow(length = unit(0.2, "cm"),ends="first"),color="black")
	}
}

# make things bigger
lollipops <- lollipops + theme(text=element_text(size=20,color="black"))

pdf(pdfFile, height = 4, width = 12)
print(lollipops)
dev.off()

}
