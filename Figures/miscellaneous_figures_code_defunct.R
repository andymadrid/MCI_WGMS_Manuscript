# Figures for LOAD / MCI / CONTROL WGMS Results

library(data.table)

################## 
################## 
# Volcano plot of DMPs
makeVolcano <- function(pvals,piGroup,title,pdf) {

	# Filter for significant/non-significant CpGs
	pvals <- as.data.frame(pvals)
	sig <- pvals[which(pvals$lfdr.from.ss<0.05 & abs(pvals[,piGroup]) >= 0.025),]
	notSig <- pvals[setdiff(rownames(pvals),rownames(sig)),]

	# Sort by hyper or hypo methylated
	sigHyper <- sig[which(sig[,piGroup]>0),]
	sigHypo <- sig[which(sig[,piGroup]<0),]

	# Take a random (25%) sample of the non-significant CpGs as to not have to plot all
	set.seed(714)
	notSig <- notSig[sample(nrow(notSig),(nrow(notSig)*0.05)),]

	# Add colors for dots
	notSig$color <- "ivory3"
	sigHyper$color <- "firebrick3"
	sigHypo$color <- "mediumblue"
	
	# Bringing it all back home
	toPlot <- rbind(sigHyper,sigHypo,notSig)

	# Plot it out…east side plot it out
	pdf(pdf)
	yl <- expression( -log [10] (lFDR))
	xl <- title
	plot(toPlot[,piGroup],-log10(toPlot$lfdr.from.ss),cex=0.7,
	     xlim=c(-0.2,0.2),col=toPlot$color,pch=16,
	     xlab=xl,
	     ylab=yl)
	abline(h=-log10(0.05),lty=2)
	abline(v=-0.025,lty=2)
	abline(v=0.025,lty=2)
	dev.off()
}
# Make volcano plots for data
df <- fread("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/MCI_CONTROL/Outputs/Summaries/pvals.bed",header=T)
makeVolcano(df,"pi.diff.mci.ctrl","Mean Methylation Difference (MCI - CU)","/media/Data/volcano.mci.control.pdf")
df <- fread("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/LOAD_CONTROL/Outputs/Summaries/pvals.bed",header=T)
makeVolcano(df,"pi.diff.load.ctrl","Mean Methylation Difference (LOAD - CU)","/media/Data/volcano.load.ctrl.pdf")
df <- fread("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/LOAD_MCI/Outputs/Summaries/pvals.bed",header=T)
makeVolcano(df,"pi.diff.load.mci","Mean Methylation Difference (LOAD - MCI)","/media/Data/volcano.laod.mci.pdf")


################## 
################## 
# Manhattan plot for DMPs
makeManhattan <- function(pvals, pdf) {
	
	# load packages needed
	library(dplyr)
	library(tidyverse)
	library(ggsci)

	# Filter for significant/non-significant CpGs
	sig <- pvals[which(pvals$lfdr.from.zz<0.05),]
	notSig <- pvals[which(pvals$lfdr.from.zz>0.05),]

	# Take a random (25%) sample of the non-significant CpGs as to not have to plot all
	set.seed(714)
	notSig <- notSig[sample(nrow(notSig),(nrow(notSig)*0.25)),]

	# Bringing it all back home
	toPlot <- rbind(sig,notSig)

	# Order chromosomes accordingly
	toPlot$chr <- factor(toPlot$chr, levels=c(paste0("chr",1:22)))
	toPlot <- toPlot[order(toPlot$chr),]

	# Calculate size of chromosomes
	toPlot2 <- toPlot %>% 
		group_by(chr) %>% 
		summarise(chr_len=max(start)) %>%

	# Calculate cumulative bases per chromosome
		mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
		select(-chr_len) %>%

	# Add data back to dataset
	left_join(toPlot, ., by=c("chr"="chr")) %>%
	arrange(chr, start) %>%
	mutate( startC=start+tot)

	# Prepare x-axis
	newAxis <- toPlot2 %>% group_by(chr) %>% summarize(center=( max(startC) + min(startC) ) / 2 )

	# Change lFDR direction depending on if hyper or hypomethylated (slow moving honey)
   	toPlot2 <- toPlot2 %>% 
		dplyr::mutate(
        new.lfdr = 
	case_when(pi.diff > 0 ~ -log10(lfdr.from.zz),
		pi.diff < 0 ~ log10(lfdr.from.zz)))

	# Prepare color palette for chromosomes
	colPal <- pal_frontiers()(9)
	colPal <- c(rep(colPal,2),colPal[1:4])

	# East side plot it out…north side plot it out
	pdf(pdf)
	man <- ggplot(toPlot2, aes(x=startC, y=new.lfdr)) +
		geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
		scale_color_manual(values = colPal) +
		scale_x_continuous( label = newAxis$chr, breaks= newAxis$center ) +
		scale_y_continuous(expand = c(0, 0) ) + 
		xlab("Chromosome") + 
		ylab("-log10(lFDR)") + 
		theme_bw() +
		theme( 
     			legend.position="none",
     			panel.border = element_blank(),
      			panel.grid.major.x = element_blank(),
      			panel.grid.minor.x = element_blank(),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(man)
	dev.off()
}
df <- fread("/media/Data/WGBS/LOAD_MCI/Results/MCI_Control/Outputs/Summaries/pvals.bed",header=T) 
makeManhattan(df,"/media/Data/manhattan.pdf")

################## 
################## 
# Make a Sankey Plot
makeSankey <- function(pvals,piGroup) {

	# Load packages
	library(networkD3)
	library(dplyr)
	library(wesanderson)
	library(ChIPseeker)
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

	# Filter for only significant CpGs
	pvals <- as.data.frame(pvals)
	sig <- pvals[which(pvals$lfdr.from.ss<0.05 & abs(pvals[,piGroup]) >= 0.025),]

	# Annotate DMPs to genes
	sig <- with(sig,GRanges(chr,IRanges(start,end)))
	peaks <- annotatePeak(sig,tssRegion=c(-3000,0),level="gene",TxDb=txdb,annoDb="org.Hs.eg.db")
	peaks <- as.data.table(peaks)

	# Change annotation to make it more streamlined
	totalNum <- nrow(peaks)
	intronPer <- 100*(length(grep("Intron",peaks$annotation))/totalNum)
	exonPer <- 100*(length(grep("Exon",peaks$annotation))/totalNum)
	promoterPer <- 100*(length(grep("Promoter",peaks$annotation))/totalNum)
	utr5Per <- 100*(length(grep("5' UTR",peaks$annotation))/totalNum)
	utr3Per <- 100*(length(grep("3' UTR",peaks$annotation))/totalNum)
	downstreamPer <- 100*(length(grep("Downstream",peaks$annotation))/totalNum)
	interPer <- 100*(length(grep("Intergenic",peaks$annotation))/totalNum)

	# Round values
	if(intronPer > 1) { intronPer <- round(intronPer)}
	if(exonPer > 1) { exonPer <- round(exonPer)}
	if(promoterPer > 1) { promoterPer <- round(promoterPer)}
	if(utr5Per > 1) { utr5Per <- round(utr5Per)}
	if(utr3Per > 1) { utr3Per <- round(utr3Per)}
	if(downstreamPer > 1) { downstreamPer <- round(downstreamPer)}
	if(interPer > 1) { interPer <- round(interPer)}

	# Adjust values if below 1%
	if(intronPer < 1 & intronPer > 0) { intronPer <- 1}
	if(exonPer < 1 & exonPer > 0) { exonPer <- 1}
	if(promoterPer < 1 & promoterPer > 0) { promoterPer <- 1}
	if(utr5Per < 1 & utr5Per > 0) { utr5Per <- 1}
	if(utr3Per < 1 & utr3Per > 0) { utr3Per <- 1}
	if(downstreamPer < 1 & downstreamPer > 0) { downstreamPer <- 1}
	if(interPer < 1 & interPer > 0) { interPer <- 1}

	# Adjust for % above 100% due to adjusting above
	x <- c("Promoter" = promoterPer, "5' UTR" = utr5Per, "Exon" = exonPer,
		"Intron" = intronPer, "3' UTR" = utr3Per, "Downstream" = downstreamPer, "Intergenic" = interPer)
	diff <- sum(x)-100
	interPer <- interPer-diff

	# set up data
	toUse <- peaks
	totalNum <- nrow(toUse)
	name1 <- paste0("Differentially Methylated Positions (DMPs) (N = ",totalNum,")")
	name2 <- paste0("Promoter (",promoterPer,"%)")
	name3 <- paste0("5' UTR (",utr5Per,"%)")
	name4 <- paste0("Exonic (",exonPer,"%)")
	name5 <- paste0("Intronic (",intronPer,"%)")
	name6 <- paste0("3' UTR (",utr3Per,"%)")
	name7 <- paste0("Downstream (",downstreamPer,"%)")
	name8 <- paste0("Intergenic (",interPer,"%)")

	# set up links
	links <- data.frame(
		source=c(name1,name1,name1,name1,name1,name1,name1),
		target=c(name2,name3,name4,name5,name6,name7,name8),

		# set up values for plot
		value=c(promoterPer,utr5Per,exonPer,intronPer,utr3Per,downstreamPer,interPer))
		nodes <- data.frame(
		name=c(as.character(links$source), 
		as.character(links$target)) %>% unique())

	# set up groups for connections
	links$group <- as.factor(c("type_a","type_b","type_c","type_d","type_e","type_f","type_g"))

	# set up groups for nodes
	nodes$group <- as.factor(c("group1"))

	# set up color scheme
	my_color <- 'd3.scaleOrdinal() .domain(["type_a","type_b","type_c","type_d","type_e","type_f","type_g","group1"]) .range(["#FF410DFF","#6EE2FFFF","#F7C530FF","#95CC5EFF","#D0DFE6FF","#F79D1EFF","#748AA6FF","black"])'
	links$IDsource <- match(links$source, nodes$name)-1 
	links$IDtarget <- match(links$target, nodes$name)-1
	sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",Value = "value", NodeID = "name",colourScale=my_color, LinkGroup="group", NodeGroup="group",fontSize=20,sinksRight=FALSE)
}

# Run locally without having to deal with server HTML issues
df <- fread("~/Desktop/pvals.mci.ctrl.bed",header=T)
makeSankey(df,"pi.diff.mci.ctrl")
df <- fread("~/Desktop/pvals.load.ctrl.bed",header=T)
makeSankey(df,"pi.diff.load.ctrl")
df <- fread("~/Desktop/pvals.load.mci.bed",header=T)
makeSankey(df,"pi.diff.load.mci")


################## 
################## 
# Gene Ontology 
runGO <- function(pvals, piGroup, pdfFile) {

	# Load package
	library(clusterProfiler)
	library(EnsDb.Hsapiens.v86)
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	library(ChIPseeker)
	EnDb <- EnsDb.Hsapiens.v86
#	seqlevelsStyle(EnDb) <- "UCSC"
#	seqlevelsStyle(EnDb) <- "NCBI"
#	seqlevelsStyle(EnDb) <- "RefSeq"
#	all.genes <- genes(EnDb)
#	autosomes <- all.genes[seqnames(all.genes) %in% 1:22]
#	autosomes$symbol[autosomes$gene_biotype == "protein_coding"]

	# Filter for only significant CpGs
	pvals <- as.data.frame(pvals)
	sig <- pvals[which(pvals$lfdr.from.ss<0.05 & abs(pvals[,piGroup]) >= 0.025),]
	sig2 <- with(sig,GRanges(chr,IRanges(start,end)))

	# Get background 
	#bg.gr <- with(pvals,GRanges(chr,IRanges(start,end)))

	# Annotate background to genes
	#peaksBG <- annotatePeak(bg.gr,tssRegion=c(-3000,0),level="gene",TxDb=txdb,annoDb="org.Hs.eg.db")
	#peaksBG <- as.data.frame(peaksBG)
	#peaksBG[which(peaksBG$annotation=="Distal Intergenic"),"SYMBOL"] <- NA
	#e <- bitr(peaksBG$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
	#genesUni <- e[,2]

	# Annotate DMPs to genes
	peaks <- annotatePeak(sig2,tssRegion=c(-3000,0),level="gene",TxDb=txdb,annoDb="org.Hs.eg.db")
	peaks <- as.data.frame(peaks)

	# Clean it up a bit
	sig <- cbind(sig,peaks$SYMBOL,peaks$distanceToTSS,peaks$annotation)
	colnames(sig) <- c(colnames(sig[,1:13]),"Symbol","distanceToTSS","annotation")
	sig[which(sig$annotation=="Distal Intergenic"),"Symbol"] <- NA
	sig$annotation <- gsub("Intron(.*)","Intron",sig$annotation)
	sig$annotation <- gsub("Exon(.*)","Exon",sig$annotation)
	sig$annotation <- gsub("Distal(.*)","Intergenic",sig$annotation)
	sig$annotation <- gsub("Promoter(.*)","Promoter",sig$annotation)
	sig$annotation <- gsub("3'(.*)","3'UTR",sig$annotation)
	sig$annotation <- gsub("5'(.*)","5'UTR",sig$annotation)
	sig$annotation <- gsub("Downstream.*)","Downstream",sig$annotation)

	# Get ENTREZID for genes
	e <- bitr(sig$Symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
	genesAll <- e[,2]

	# run GOs
	#egoAll <- enrichGO(gene=genesAll,universe=genesUni,ont="BP",
	#	readable=T,OrgDb="org.Hs.eg.db",pvalueCutoff=0.8,qvalueCutoff=1)
	egoAll <- enrichGO(gene=genesAll,ont="BP",
		readable=T,OrgDb="org.Hs.eg.db",pvalueCutoff=0.8,qvalueCutoff=1)
	egoAll <- simplify(egoAll)

	# load packages for plotting
	library(forcats)
	library(gplots)
	library(ggplot2)
	library(wesanderson)
	colPal <- wes_palette("Zissou1",25,"continuous")
	egoAll <- as.data.frame(egoAll)
	egoAll$Description <- factor(egoAll$Description,levels=egoAll$Description)
	pdf(pdfFile,height=8,width=12)
	p1 <- ggplot(egoAll[1:25,],aes(x=-log10(pvalue),y=fct_rev(Description),fill=Description)) + 
		geom_bar(stat="identity") + 
		scale_fill_manual(values=colPal) + 
		theme_classic() + 
		theme(legend.position="none") + 
		ggtitle("Differentially Methylated Genes") + 
		theme(axis.title.y=element_blank(), text=element_text(size=20))
	print(p1)
	dev.off()

	# get GO results out
	return(egoAll)
}
df <- fread("~/Desktop/pvals.mci.ctrl.bed",header=T)
egoResults <- runGO(df,"pi.diff.mci.ctrl","~/Desktop/go.mci.ctrl.pdf")
df <- fread("~/Desktop/pvals.load.ctrl.bed",header=T)
egoResults <- runGO(df,"pi.diff.load.ctrl","~/Desktop/go.load.ctrl.pdf")
df <- fread("~/Desktop/pvals.load.mci.bed",header=T)
egoResults <- runGO(df,"pi.diff.load.mci","~/Desktop/go.load.mci.pdf")


################## 
################## 
# Density plots
makeDense <- function(pdf) {

	# load packages
	library(dmrseq)
	library(ggplot2)
#	source("https://github.com/andymadrid/One-Off-Scripts/blob/main/plotEmpiricalDistribution2.r")

	# read in phenotypic data
	targets <- read.csv("/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv",header=T)

	# change directory to where beds are
	setwd("/media/Data/WGBS/LOAD_MCI/MCovRaw/BismarkFormatted/")
	infile <- list.files(pattern="bed$")

	# filter out LOAD data
	infile <- infile[which(targets$diagnostic_group!="LOAD")]
	targets <- targets[which(targets$diagnostic_group!="LOAD"),]

	# read in bed files
	bs <- read.bismark(files = infile,rmZeroCov=TRUE,strandCollapse=T,verbose=T)

	# filter out CpGs with zero coverage in all samples for the sake of plotting
	loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0) == 0)
	bs <- bs[loci.idx,]

	# take a random sample (1 million) of CpGs for easier plotting - considering large sample size
	set.seed(714)
	bs <- bs[sample(nrow(bs),1e6),]
		
	# Add metadata
	pData(bs)$Condition <- targets$diagnostic_group

	# plot it out
	pdf(pdf)
	plotEmpiricalDistribution2(bs, testCovariate = "Condition", bySample = FALSE, type = "M", colPal=c("mediumblue","firebrick3"))
	plotEmpiricalDistribution2(bs, testCovariate = "Condition", bySample = FALSE, type = "Cov", colPal=c("mediumblue","firebrick3"))
	dev.off()
}
makeDense("/media/Data/densityPlots.pdf")


################## 
################## 
# Heatmap of DMPs (using pis)

plotHeat <- function(pvals, pis.dir, pdf) {
	
	# load packages
	cat("Loading in packages now . . .\n")
	suppressPackageStartupMessages({
	library(gplots)
	library(RColorBrewer)
	library(dmrseq)
	library(viridis)
	library(data.table)
	library(GenomicRanges)
	library(dplyr)
	})

	# read in estimated methylation levels
	cat("Reading in pi data . . . this may take awhile . . . go grab a beer or something . . .\n\n")
	setwd(pis.dir)
	pi.files <- list.files(pattern = "bed$")
	pi.df <- lapply(pi.files, fread)
	pi.df <- as.data.frame(do.call(rbind,pi.df))

	# read in phenotypic data
	targets <- read.csv("/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv",header=T)
	targets <- targets[match(colnames(pi.df), targets$sample_id),]
	targets <- targets[-c(1:3),]
	targets <- targets %>%
		mutate(color = 
		case_when(diagnostic_group == "MCI" ~ "darkgoldenrod2",
		diagnostic_group == "CONTROL" ~ "dodgerblue2"))

	# filter for only significant CpGs
	cat("Filtering for DMPs . . .\n")
	dmps <- pvals[which(pvals$lfdr.from.ss < 0.05),]

	# get chromosome and position information
	cat("Creating granges objects . . .\n")
	meth.gr <- with(pi.df,GRanges(chr,IRanges(start)))
	dmps.gr	<- with(dmps,GRanges(chr,IRanges(start)))

	# Filter pis for only DMPs
	cat("Finding overlaps and filtering objects . . .\n")
	overlaps <- as.data.frame(findOverlaps(meth.gr, dmps.gr))
	pi.sub <- pi.df[overlaps[,1], ]
	pi.sub <- pi.sub[,-c(1:3)]

	# make color palette
	cat("Look at all of the pretty colors!!\n")
	pal <- colorpanel(100,"dodgerblue3","goldenrod1","firebrick3")

	# Plot it out . . . south side plot it out
	cat("Plotting . . .\n")
	pdf(pdf)
	heatmap.2(as.matrix(pi.sub),col=pal,trace="none",
		dendrogram="column",labRow=F,labCol=F,key.title=NULL,
		key.xlab="Methylation Level",density.info="none",
		ColSideColors=targets$color)
	dev.off()
}
df <- fread("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/MCI_CONTROL/Outputs/Summaries/pvals.bed",header=T) 
plotHeat(df, "/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/MCI_CONTROL/Outputs/Pis/", "/media/Data/heatmap.pdf")



################## 
################## 
# Violin plots

# load in some packages
library(data.table)
library(ggplot2)
library(GenomicRanges)

# read in differential results
df <- fread("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/CONTROL_MCI_LOAD/Outputs/Summaries/pvals.bed",header=T)

# filter to significant CpGs
#sig <- df[which(df$lfdr.from.ss<0.05),] 
sig <- df[which(df$lfdr.from.ss<0.05 & ((df$pi.diff.mci.ctrl > 0.0125 & df$pi.diff.load.mci > 0.0125) | (df$pi.diff.mci.ctrl < -0.0125 & df$pi.diff.load.mci < -0.0125))),] 
sig <- sig[order(sig$p.from.DSS),]
setwd("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/CONTROL_MCI_LOAD/Outputs/Pis/")

# load in estimated methylation values data
pi.files <- list.files(pattern = "bed$")
pi.df <- lapply(pi.files, fread)
pi.df <- as.data.frame(do.call(rbind,pi.df))
meth.gr <- with(pi.df,GRanges(chr,IRanges(start)))
dmps.gr	<- with(sig,GRanges(chr,IRanges(start)))
overlaps <- as.data.frame(findOverlaps(dmps.gr, meth.gr))
pi.sub <- pi.df[overlaps[,2], ]
cpg.info <- paste0(pi.sub$chr,":",pi.sub$start+1)
pi.sub <- pi.sub[,-c(1:3)]
pi.sub <- as.data.frame(t(pi.sub))
pi.sub$groups <- gsub("LOAD", "AD", pi.sub$groups)

# load in phenotypic data
targets <- read.csv("/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv",header=T)
#removeIDs <- read.table("/media/Data/WGBS/LOAD_MCI/removeSamples.txt",header=T)
#targets <- subset(targets,!(study_id %in% removeIDs$sampleID))

# Bringing it all back home
pi.sub$groups <- targets$diagnostic_group

# Create dataset to plot
i <- 1 # first hypo DMP
#i <- 5 # first hyper DMP
ii <- ncol(pi.sub)
mainLab <- cpg.info[i]
x <- pi.sub[,c(i,ii)]
colnames(x) <- c("pis", "groups")
x$groups <- factor(x$groups, levels = c("CONTROL", "MCI", "AD"))

# Plot it out
pdf("/media/Data/vln.pdf")
p <- ggplot(x, aes(x=groups, y=pis*100, fill=groups)) +
	geom_violin(trim=F) +
	scale_fill_manual(values=c("firebrick3","darkgoldenrod3","mediumblue")) +
	#geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25) +
	theme_bw() +
	xlab("") +
	ylab("Estimated Methylation Level (%)") + labs(title=mainLab) +
	theme(legend.position="none") +
	theme(text=element_text(size=20,color="black"))
print(p)
dev.off()




################## 
################## 
# Machine learning of DMPs
machLearn <- function(pvals, pis.dir, pdf) {

	# load packages
	suppressPackageStartupMessages({
	library(DMRichR)
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	library(org.Hs.eg.db)
	library(dplyr)
	library(gplots)
	library(RColorBrewer)
	library(dmrseq)
	library(viridis)
	library(data.table)
	library(GenomicRanges)
	})

	# read in estimated methylation levels
	cat("Reading in pi data . . . this may take awhile . . . go grab a beer or something . . .\n\n")
	setwd(pis.dir)
	pi.files <- list.files(pattern = "bed$")
	pi.df <- lapply(pi.files, fread)
	pi.df <- as.data.frame(do.call(rbind,pi.df))

	# filter for only significant CpGs
	cat("Filtering for DMPs . . .\n")
	dmps <- pvals[which(pvals$lfdr.from.zz < 0.05),]

	# get chromosome and position information
	cat("Creating granges objects . . .\n")
	meth.gr <- with(pi.df,GRanges(chr,IRanges(start)))
	dmps.gr	<- with(dmps,GRanges(chr,IRanges(start)))

	# Filter pis for only DMPs
	cat("Finding overlaps and filtering objects . . .\n")
	overlaps <- as.data.frame(findOverlaps(meth.gr, dmps.gr))
	pi.sub <- pi.df[overlaps[,1], ]
	pi.sub <- pi.sub[,-c(1:3)]
	pi.sub <- t(as.matrix(pi.sub))
	colnames(pi.sub) <- cbind(as.data.frame(seqnames(dmps.gr)), ranges(dmps.gr)) %>%
		dplyr::as_tibble() %>%
		dplyr::select(value,start) %>%
		tidyr::unite("bed", c("value","start"), sep = ".") %>%
		as.matrix()

	# read in phenotypic data
	targets <- read.csv("/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv",header=T)
	targets <- targets[match(rownames(pi.sub), targets$sample_id),]
	targets <- targets %>%
		mutate(color = 
		case_when(diagnostic_group == "MCI" ~ "red",
		diagnostic_group == "CONTROL" ~ "blue"))

	# Add group information to data
	pi.sub <- pi.sub %>%
		dplyr::as_tibble() %>%
		tibble::add_column(groups = targets$diagnostic_group) %>%
		tibble::add_column(sampleID = rownames(pi.sub))

	# Helper function to split DMPs to chr and start and add to tibble
	splitDmrs <- function(ranking) {
		DMRsplit <- ranking$DMP %>%
			strsplit(., split = "[.]") %>%
			as.data.frame() %>%
			t() %>%
			magrittr::set_colnames(c("chr", "start")) %>%
			dplyr::as_tibble()
    
		ranking <- ranking %>%
			tibble::add_column(chr = DMRsplit$chr, .after = 2) %>%
			tibble::add_column(start = DMRsplit$start, .after = 3)

		return(ranking)
	}
		
	# Run support vector machine learning algorithm
	cat("Running support vector machine RFE algorithm now . . .\n")
	dataMatrix <- pi.sub %>%
		dplyr::select(-c(groups, sampleID)) %>% 
		as.matrix()
	set.seed(714)
	svmTrained <- sigFeature::sigFeature(dataMatrix, pi.sub$groups)

	svmRanked <- tibble::tibble(Rank = 1:length(svmTrained),
		DMP = colnames(dataMatrix[, svmTrained]))
	
	svmRanked <- svmRanked %>% splitDmrs()
		
	# Run random forest algorithm
	cat("Running random forest algorithm now . . .\n")
	set.seed(714)
	rfTrained <- Boruta::Boruta(factor(groups) ~ ., data = pi.sub %>% dplyr::select(-sampleID), doTrace = 0)
	rfStats <- Boruta::attStats(rfTrained)
	rfRanking <- tibble::tibble(DMP = rownames(rfStats),
		meanImp = rfStats$meanImp,
		decision = rfStats$decision) %>%
	dplyr::arrange(dplyr::desc(meanImp)) %>%
	tibble::add_column(Rank = 1:nrow(rfStats), .before = 1)
	rfRanking <- rfRanking %>% splitDmrs()

	# Find common DMPs among the top 1% of ranked DMPs by the two algorithms
	cat("Finding common DMPs among the top 1% of the algorithms . . .\n")
	nPredictors <- ncol(pi.sub) - 2
	nTopPercent <- ceiling(1 * .01 * nPredictors)
	commonDMPs <- intersect(rfRanking$DMP[1:nTopPercent], svmRanked$DMP[1:nTopPercent])
	commonDMPsRfRank <- which(rfRanking$DMP %in% commonDMPs)
	commonDMPsSvmRank <- which(svmRanked$DMP %in% commonDMPs)
	case <- nTopPercent
	#commonDMPs <- list(DMP = commonDMPs, rfRank = commonDMPsRfRank , svmRank = commonDMPsSvmRank, case = case)

	# Filter pis to only the most predictive ones
	commonDMPs <- intersect(rfRanking$DMP[1:nTopPercent], svmRanked$DMP[1:nTopPercent])
	commonDMPs <- as.data.frame(commonDMPs)
	commonDMPs <- as.data.frame(do.call(rbind,commonDMPs$commonDMPs %>% strsplit(., split = "[.]")))
	colnames(commonDMPs) <- c("chr","start")
	predictiveDMPs <- c()
	predictivePis <- c()
	pi.sub <- pi.df[overlaps[,1], ]
	for (i in 1:nrow(commonDMPs)) {
	predictiveDMPs <- rbind(predictiveDMPs, dmps[dmps$chr==commonDMPs[i,"chr"] & dmps$start==commonDMPs[i,"start"],])
	predictivePis <- rbind(predictivePis, pi.sub[pi.sub$chr==commonDMPs[i,"chr"] & pi.sub$start==commonDMPs[i,"start"],])
	}
	predictivePis <- predictivePis[,-c(1:3)]

	# make color palette
	cat("Look at all of the pretty colors!!\n")
	pal <- colorpanel(100,"dodgerblue3","goldenrod1","firebrick3")

	# Plot a heat map of the most predictive DMPs
	cat("Making a heat map of the most predictive DMPs . . .\n")
	pdf(pdf)
	heatmap.2(as.matrix(predictivePis),col=pal,trace="none",
		dendrogram="column",labRow=F,labCol=F,key.title=NULL,
		key.xlab="Methylation Level",density.info="none",
		ColSideColors=targets$color)
	dev.off()

}

df <- fread("/media/Data/WGBS/LOAD_MCI/Results/MCI_Control/Outputs/Summaries/pvals.bed",header=T) 
machLearn(df,"/media/Data/WGBS/LOAD_MCI/Results/MCI_Control/Outputs/Pis/","/media/Data/heatmap.predictiveDMPs.pdf")




################## 
################## 
# Machine learning of DMPs (CONTROL / MCI / LOAD)
# Using training (80%) and validation sets (20%)
machLearn <- function(pvals, pis.dir, pdfFile, pdfFile2) {

	# load packages
	suppressPackageStartupMessages({
	library(DMRichR)
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	library(org.Hs.eg.db)
	library(dplyr)
	library(gplots)
	library(RColorBrewer)
	library(dmrseq)
	library(viridis)
	library(data.table)
	library(GenomicRanges)
	})

	# read in estimated methylation levels
	cat("Reading in pi data . . . this may take awhile . . . go grab a beer or something . . .\n\n")
	setwd(pis.dir)
	pi.files <- list.files(pattern = "bed$")
	pi.df <- lapply(pi.files, fread)
	pi.df <- as.data.frame(do.call(rbind,pi.df))

	# filter for only significant CpGs
	cat("Filtering for DMPs . . .\n")
	dmps <- pvals[which(pvals$lfdr.from.ss < 0.05),]

	# get chromosome and position information
	cat("Creating granges objects . . .\n")
	meth.gr <- with(pi.df,GRanges(chr,IRanges(start)))
	dmps.gr	<- with(dmps,GRanges(chr,IRanges(start)))

	# Filter pis for only DMPs
	cat("Finding overlaps and filtering objects . . .\n")
	overlaps <- as.data.frame(findOverlaps(meth.gr, dmps.gr))
	pi.sub <- pi.df[overlaps[,1], ]
	pi.sub <- pi.sub[,-c(1:3)]
	pi.sub <- t(as.matrix(pi.sub))
	colnames(pi.sub) <- cbind(as.data.frame(seqnames(dmps.gr)), ranges(dmps.gr)) %>%
		dplyr::as_tibble() %>%
		dplyr::select(value,start) %>%
		tidyr::unite("bed", c("value","start"), sep = ".") %>%
		as.matrix()

	# read in phenotypic data
	targets <- read.csv("/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv",header=T)
	targets <- targets[match(rownames(pi.sub), targets$sample_id),]
	targets <- targets %>%
		mutate(color = 
		case_when(diagnostic_group == "MCI" ~ "darkgoldenrod2",
		diagnostic_group == "LOAD" ~ "firebrick2",
		diagnostic_group == "CONTROL" ~ "dodgerblue2"))
	set.seed(714)
	targets.train <- targets[sample(nrow(targets),(nrow(targets)*0.8)),]
	valSamples <- setdiff(targets$sample_id,targets.train$sample_id)
	targets.validation <- c()
	for (i in 1:length(valSamples)) {
		targets.validation <- rbind(targets.validation,targets[which(targets$sample_id==valSamples[i]),])
	}
	pi.sub.train <- pi.sub[match(targets.train$sample_id, rownames(pi.sub)),]
	pi.sub.validation  <- pi.sub[match(targets.validation$sample_id, rownames(pi.sub)),]

	# Add group information to data
	pi.sub.train <- pi.sub.train %>%
		dplyr::as_tibble() %>%
		tibble::add_column(groups = targets.train$diagnostic_group) %>%
		tibble::add_column(sampleID = rownames(pi.sub.train))
	pi.sub.validation <- pi.sub.validation %>%
		dplyr::as_tibble() %>%
		tibble::add_column(groups = targets.validation$diagnostic_group) %>%
		tibble::add_column(sampleID = rownames(pi.sub.validation))

	# Helper function to split DMPs to chr and start and add to tibble
	splitDmrs <- function(ranking) {
		DMRsplit <- ranking$DMP %>%
			strsplit(., split = "[.]") %>%
			as.data.frame() %>%
			t() %>%
			magrittr::set_colnames(c("chr", "start")) %>%
			dplyr::as_tibble()
    
		ranking <- ranking %>%
			tibble::add_column(chr = DMRsplit$chr, .after = 2) %>%
			tibble::add_column(start = DMRsplit$start, .after = 3)

		return(ranking)
	}
		
	# Run support vector machine learning algorithm
	cat("Running support vector machine RFE algorithm on training set now . . .\n")
	dataMatrix <- pi.sub.train %>%
		dplyr::select(-c(groups, sampleID)) %>% 
		as.matrix()
	set.seed(714)
	svmTrained <- sigFeature::sigFeature(dataMatrix, pi.sub.train$groups)

	svmRanked <- tibble::tibble(Rank = 1:length(svmTrained),
		DMP = colnames(dataMatrix[, svmTrained]))
	
	svmRanked <- svmRanked %>% splitDmrs()
		
	# Run random forest algorithm
	cat("Running random forest algorithm on training set now . . .\n")
	set.seed(714)
	rfTrained <- Boruta::Boruta(factor(groups) ~ ., data = pi.sub.train %>% dplyr::select(-sampleID), doTrace = 0)
	rfStats <- Boruta::attStats(rfTrained)
	rfRanking <- tibble::tibble(DMP = rownames(rfStats),
		meanImp = rfStats$meanImp,
		decision = rfStats$decision) %>%
	dplyr::arrange(dplyr::desc(meanImp)) %>%
	tibble::add_column(Rank = 1:nrow(rfStats), .before = 1)
	rfRanking <- rfRanking %>% splitDmrs()

	# Find common DMPs among the top 1% of ranked DMPs by the two algorithms
	cat("Finding common DMPs among the top 1% of the algorithms . . .\n")
	nPredictors <- ncol(pi.sub.train) - 2
	nTopPercent <- ceiling(1 * .01 * nPredictors)
	commonDMPs <- intersect(rfRanking$DMP[1:nTopPercent], svmRanked$DMP[1:nTopPercent])
	commonDMPsRfRank <- which(rfRanking$DMP %in% commonDMPs)
	commonDMPsSvmRank <- which(svmRanked$DMP %in% commonDMPs)
	case <- nTopPercent
	#commonDMPs <- list(DMP = commonDMPs, rfRank = commonDMPsRfRank , svmRank = commonDMPsSvmRank, case = case)

	# Filter pis to only the most predictive ones
	commonDMPs <- intersect(rfRanking$DMP[1:nTopPercent], svmRanked$DMP[1:nTopPercent])
	commonDMPs <- as.data.frame(commonDMPs)
	commonDMPs <- as.data.frame(do.call(rbind,commonDMPs$commonDMPs %>% strsplit(., split = "[.]")))
	colnames(commonDMPs) <- c("chr","start")
	predictiveDMPs.train <- c()
	predictivePis.train <- c()
	pi.sub.train2 <- pi.df[overlaps[,1], ]
	x <- pi.sub.train2[,c("chr", "start", "end")]
	pi.sub.train2 <- pi.sub.train2[,match(targets.train$sample_id, colnames(pi.sub.train2))]
	pi.sub.train2 <- cbind(x,pi.sub.train2)
	for (i in 1:nrow(commonDMPs)) {
	predictiveDMPs.train <- rbind(predictiveDMPs.train,
		dmps[dmps$chr==commonDMPs[i,"chr"] & dmps$start==commonDMPs[i,"start"],])
	predictivePis.train <- rbind(predictivePis.train,
		pi.sub.train2[pi.sub.train2$chr==commonDMPs[i,"chr"] & pi.sub.train2$start==commonDMPs[i,"start"],])
	}
	predictivePis.train <- predictivePis.train[,-c(1:3)]
	predictivePis.train <- data.matrix(predictivePis.train)

	predictiveDMPs.validation <- c()
	predictivePis.validation <- c()
	pi.sub.validation2 <- pi.df[overlaps[,1], ]
	x <- pi.sub.validation2[,c("chr", "start", "end")]
	pi.sub.validation2 <- pi.sub.validation2[,match(targets.validation$sample_id, colnames(pi.sub.validation2))]
	pi.sub.validation2 <- cbind(x,pi.sub.validation2)
	for (i in 1:nrow(commonDMPs)) {
	predictiveDMPs.validation <- rbind(predictiveDMPs.validation,
		dmps[dmps$chr==commonDMPs[i,"chr"] & dmps$start==commonDMPs[i,"start"],])
	predictivePis.validation <- rbind(predictivePis.validation,
		pi.sub.validation2[pi.sub.validation2$chr==commonDMPs[i,"chr"] & pi.sub.validation2$start==commonDMPs[i,"start"],])
	}
	predictivePis.validation <- predictivePis.validation[,-c(1:3)]
	predictivePis.validation <- data.matrix(predictivePis.validation)

	# make color palette
	cat("Look at all of the pretty colors!!\n")
	pal <- colorpanel(100,"dodgerblue3","goldenrod1","firebrick3")

	# Plot a heat map of the most predictive DMPs
	cat("Making a heat map of the most predictive DMPs . . .\n")
	pdf(pdfFile)
	heatmap.2(as.matrix(predictivePis.train),col=pal,trace="none",
		dendrogram="column",labRow=F,labCol=F,key.title=NULL,
		key.xlab="Methylation Level",density.info="none",
		ColSideColors=targets.train$color)
	dev.off()

	pdf(pdfFile2)
	heatmap.2(as.matrix(predictivePis.validation),col=pal,trace="none",
		dendrogram="column",labRow=F,labCol=F,key.title=NULL,
		key.xlab="Methylation Level",density.info="none",
		ColSideColors=targets.validation$color)
	dev.off()

}

df <- fread("/media/Data/WGBS/LOAD_MCI/Results/MCI_Control/Outputs/Summaries/pvals.bed",header=T) 
machLearn(df,"/media/Data/WGBS/LOAD_MCI/Results/MCI_Control/Outputs/Pis/","/media/Data/heatmap.train.pdf","/media/Data/heatmap.validation.pdf")


########### 
########### 
## Find coverage/% imputed for predictive DMPs
## Slow process…could be a better way, but lazy to figure it out

# load packages 
suppressPackageStartupMessages({
    library(data.table)
    library(argparse)
    library(DSS)
    library(magrittr)
    library(bsseq)
    library(dplyr)
})


# loop by chromosome
# get chromosome list
chrs <- list.files("/media/Data/WGBS/LOAD_MCI/MCovRaw/BismarkFormatted/chromosomes")
chrs <- chrs[1:22]
loopByChromosome0 <- function(chr) {
cumM <- c()
cumCov <- c()
for (i in 1:length(chr)) {
cat(chr[[i]])
cat("\n")

parser <- ArgumentParser()
parser$add_argument("--idir", default= "/media/Data/WGBS/LOAD_MCI/MCovRaw/", help='Directory to run DSS on')
parser$add_argument("--master_file", default="/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv")
#parser$add_argument("--chr", default= "chr22", help='Chromosome to run on')
parser$add_argument("--chr", default= chr[[i]], help='Chromosome to run on')
parser$add_argument("--save_pivots", action="store_true", help= "Save pivoted M and Cov matrices? Usually used when you want to check.")
args <- parser$parse_args()

master.df <- read.csv(args$master_file) %>%
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::select(-starts_with("PC")) %>%
    dplyr::mutate(
        diagnostic_group_coded = 
	case_when(diagnostic_group == "CONTROL" ~ 0,
	diagnostic_group == "MCI" ~ 1,
	diagnostic_group == "LOAD" ~ 2))

# Check that the samples in the directory have phenotype and vice versa
idir.samples <- 
    unique(stringr::str_split_fixed(list.files(args$idir), "\\.", 3)[ ,2])

# All samples that are in the directory and sampleshseet
valid.samples <- intersect(master.df$sample_id, idir.samples)
load.samples <- intersect(
    master.df$sample_id[master.df$diagnostic_group == "LOAD"], 
    idir.samples)

mci.samples <- intersect(
    master.df$sample_id[master.df$diagnostic_group == "MCI"], 
    idir.samples)

ctrl.samples <- intersect(
    master.df$sample_id[master.df$diagnostic_group == "CONTROL"], 
    idir.samples
    )


read_wrapper <- function(s){
    #Creates input name like ../Data/chr22.100.bed
    # and reads it
    dt <- fread(file.path(args$idir, paste0(args$chr, ".", s ,".bed")))
    dt$sample <- as.numeric(s) # numeric quiets warning
    dt
}

data <- do.call(rbind, lapply(X=valid.samples, FUN=read_wrapper))
colnames(data) <- c("chrom","chromStart","chromEnd","strand","methylated","coverage","sample")

pivot_me <- function(data, value) {
    # Gets into Sample \times Position matrix of 
    # M or Cov
    # "value" is "methylated" or "coverage"
    keepcols <- c("chromStart", "sample", value)

    data %>% 
        dplyr::select(all_of(keepcols)) %>%
        tidyr::pivot_wider(values_from = value, names_from = sample) %>%
        tibble::column_to_rownames("chromStart")

 }

M <- pivot_me(data, "methylated")
rownames(M) <- paste0(chr[[i]],".",rownames(M))
Cov <- pivot_me(data, "coverage")
rownames(Cov) <- paste0(chr[[i]],".",rownames(Cov))
cumM <- rbind(cumM,M)
cumCov <- rbind(cumCov,Cov)
}
fwrite(cumM,file="/media/Data/cummulativeM.txt",sep="\t",row.names=T)
fwrite(cumCov,file="/media/Data/cummulativeCov.txt",sep="\t",row.names=T)
}

# load predictive DMPs
load("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/predictiveDMPs.rdata")
commonDMPs$look <- paste0(commonDMPs$chr,".",commonDMPs$start)

# load cumulative coverage file
cumCov <- fread("/media/Data/cummulativeCov.txt",header=T)
which(cumCov$V1==commonDMPs$look)


############### 
############### 
### Predicted Age
############### 
############### 

setwd("/media/Data/WGBS/LOAD_MCI/MCovRaw/BismarkFormatted/")
library(wateRmelon)
library(dmrseq)
library(DunedinPACE)
library(devtools)
library(data.table)
library(dplyr)
library(DunedinPoAm38)
infile <- list.files(pattern="bed$")
bs <- read.bismark(files = infile,rmZeroCov=F,strandCollapse=T,verbose=T)
M <- getCoverage(bs,type="M")
Cov <- getCoverage(bs,type="Cov")
chr <- as.data.frame(seqnames(bs))
colnames(chr) <- "chr"
pos <- as.data.frame(start(bs))
colnames(pos) <- "pos"
P <- as.data.frame(M/Cov)
P$chr <- chr
P$pos <- pos
P$pos <- as.numeric(unlist(P$pos))
P$chr <- as.character(unlist(P$chr))
beta.gr <- with(P,GRanges(chr,IRanges(pos,pos)))
epic <- read.table("/media/Data/AMC_WGBS/Analysis_C_v_PC/hglft_genome_1244c_b431c0.bed",header=F)
epic.gr <- with(epic,GRanges(V1,IRanges(V2)))
overlaps <- as.data.frame(findOverlaps(epic.gr,beta.gr))
epic.subset <- epic[overlaps$queryHits,]
beta.subset <- P[overlaps$subjectHits,]
rownames(beta.subset) <- epic.subset$V5
beta.subset <- as.data.frame(beta.subset[,1:length(infile)])
beta.seq <- as.matrix(beta.subset)
beta.seq[is.na(beta.seq)] <- 0
predAge.seq <- agep(beta.seq, method='all')
dune.seq <- PACEProjector(beta.seq)
targets <- read.csv("/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv",header=T)
predAge.seq$chronAge <- targets$age_at_visit
predAge.seq$diff.horvath <- predAge.seq$horvath.horvath.age - predAge.seq$chronAge
predAge.seq$diff.hannum <- predAge.seq$hannum.hannum.age - predAge.seq$chronAge
predAge.seq$diff.phenoage <- predAge.seq$phenoage.phenoage.age - predAge.seq$chronAge
predAge.seq$diff.skinblood <- predAge.seq$skinblood.skinblood.age - predAge.seq$chronAge
predAge.seq$diff.li <- predAge.seq$lin.lin.age - predAge.seq$chronAge
predAge.seq$group <- targets$diagnostic_group
predAge.seq <- predAge.seq %>%
    dplyr::mutate(
        color = 
	case_when(group == "CONTROL" ~ "#E41A1C",
	group == "MCI" ~ "#4DAF4A",
	group == "LOAD" ~ "#377EB8"))
save(predAge.seq,dune.seq,beta.seq,file="/media/Data/WGBS/LOAD_MCI/predictedAges.rdata")

# Filter to just the LOAD and CONTROL samples
predAge.seq <- predAge.seq[which(predAge.seq$group!="MCI"),]

# Plot correlation between predicted and actual age
pdf("/media/Data/cor.pdf")
plot(predAge.seq$chronAge,predAge.seq$horvath.horvath.age,pch=16,col=predAge.seq$color,xlab="Chronological Age (Years)",ylab="Estimated Age (Horvath)")
legend("topleft", legend=c("Control", "LOAD"), pch = 16, col = c("#E41A1C", "#377EB8"))
abline(0,1,col="black")
abline(lm(predAge.seq$horvath.horvath.age ~ predAge.seq$chronAge),col="black",lty=2)
dev.off()
pdf("/media/Data/bar.together.pdf")
plot(factor(predAge.seq$group),predAge.seq$diff.horvath,col=c("#E41A1C","#377EB8"),xlab="",ylab="Difference in Age (Predicted - Chronological)")
dev.off()

# Plot differences versus group
pdf("/media/Data/boxplot.pdf")
x <- ggplot(data=predAge.seq,aes(x=group,y=diff.horvath)) +  geom_boxplot(fill=c("#E41A1C","#377EB8"),outlier.shape = NA) + geom_jitter(color="black", size=0.4, alpha=0.9) + theme_classic() + theme(legend.position="none",text = element_text(size=20,face="bold")) + xlab("") + ylab("Difference in Age (Predicted - Chronological)")
x <- x + geom_signif(comparisons = list(c("CONTROL", "LOAD")),test="t.test")
print(x)
dev.off()

#pdf("/media/Data/cor.pdf")
#plot(predAge.seq$chronAge,predAge.seq$skinblood.skinblood.age,pch=16,col=predAge.seq$color,xlab="Chronological Age (Years)",ylab="Estimated Age (Skin-Blood)")
#legend("topleft", legend=c("Control", "MCI", "LOAD"), pch = 16, col = c("#E41A1C", "#377EB8", "#4DAF4A"))
#abline(0,1,col="black")
#abline(lm(predAge.seq$skinblood.skinblood.age ~ predAge.seq$chronAge),col="black",lty=2)
#dev.off()
