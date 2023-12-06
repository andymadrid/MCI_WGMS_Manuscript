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

# get nature genes (GWAS)
natgen.symbols <- get_natgen_genes("/media/Data/WGBS/LOAD_MCI/Results/adRiskGenes.txt",exclude_APP = T,exclude_IGH = T)

# read in data
pvals.data <- get_pvals_data("/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/MCI_CONTROL/Outputs/Summaries/pvals.bed")

# wrangle DMPs
pvals.gr <- to_granges(pvals.data)
dmps.data <- subset_dmps(pvals.data)
dmps.gr <- to_granges(dmps.data)

# annotate DMPs
dmps.genic.df <- annotate_loci_to_genic_parts(dmps.gr)
dmps.cpg.df <- cpg_annotation_routine(dmps.gr)
plot_cpg_pi_chart(dmps.cpg.df,"/media/Data/piMCI.pdf")
screenshot_sankey(plot_sankey(dmps.genic.df, nrow(dmps.data)),"/media/Data/dmps-genic-sankey.png")

# check overlap with nature GWAS genes
genes.expanded.25kb <- get_gene_bodies(upstream = 25000, downstream = 25000, autosomes_only = T,  protein_coding_only = T)
NG.genes.expanded.25kb <- filter_by_gene_symbols(genes.expanded.25kb, natgen.symbols)
NG.genes.with.dmp <- make_df_from_two_overlapping_granges(dmps.gr, NG.genes.expanded.25kb)
NG.tally.df <- tally_dmps_in_genes(NG.genes.with.dmp)

# annotate DMPs to gene bodies
gene.bodies <- get_gene_bodies(upstream=3000,downstream=200,autosomes_only=T,protein_coding_only=T)
dmps.in.gene.stats <- tally_dmps_in_out_genes(dmps.gr, gene.bodies)
gene.bodies.with.dmp <- make_df_from_two_overlapping_granges(dmps.gr, gene.bodies)
gene.bodies.tally.df <- tally_dmps_in_genes(gene.bodies.with.dmp)
N.DMPs.not.in.genes <- get_N_not_in_set(dmps.gr, gene.bodies)

# look at promoters only
promoters <- get_protein_coding_promoters(upstream = 5000, downstream = 200)
gene.body.enrichment <- harmonic_pvalue_routine(pvals.gr, gene.bodies, ALPHA)
gene.body.enrichment.hyper.hypo <- harmonic_pvalue_routine_hyper_hypo(pvals.gr, gene.bodies, ALPHA)
promoter.enrichment <- harmonic_pvalue_routine(pvals.gr, promoters, ALPHA)
promoter.enrichment.hyper.hypo <- harmonic_pvalue_routine_hyper_hypo(pvals.gr, promoters, ALPHA)

# clean’r up
#dm.genes.df <- dplyr::filter(gene.body.enrichment, lfdr < DMGENE.ALPHA,N.CpGs > 0)
dm.genes.df <- dplyr::filter(gene.body.enrichment, lfdr < DMGENE.ALPHA,N.DMPs > 0)
dm.genes.df.hyper <- dplyr::filter(gene.body.enrichment.hyper.hypo$hyper, lfdr < DMGENE.ALPHA,N.DMPs > 0)
dm.genes.df.hypo <- dplyr::filter(gene.body.enrichment.hyper.hypo$hypo, lfdr < DMGENE.ALPHA,N.DMPs > 0)

# gene ontology
dm.genes.go.df <- symbols_to_gene_ontology_routine(dm.genes.df$gene_name)
dm.genes.go.df.hyper <- symbols_to_gene_ontology_routine(dm.genes.df.hyper$gene_name)
dm.genes.go.df.hypo <- symbols_to_gene_ontology_routine(dm.genes.df.hypo$gene_name)
plot_go_barchart(dm.genes.go.df, n = 25,"/media/Data/go.DMGenes.pdf")
plot_go_barchart(dm.genes.go.df.hyper, n = 25,"/media/Data/go.DMGenes.hyper.pdf")
plot_go_barchart(dm.genes.go.df.hypo, n = 25,"/media/Data/go.DMGenes.hypo.pdf")
save(dm.genes.df,dm.genes.df.hyper,dm.genes.df.hypo,file="/media/Data/goData.rdata")

# PCHiC Analysis
chain.19to38 <- download_chain_from_ucsc("/media/Data/WGBS/LOAD_MCI/","https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")
chain.38to19 <- download_chain_from_ucsc("/media/Data/WGBS/LOAD_MCI/","https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz")
interactions.hg19 <- clean_interactions_data("/media/Data/WGBS/LOAD_MCI/PCHiC_peak_matrix_cutoff5.tsv")
interactions.hg38 <- lift_promoter_capture_data_to_hg38(interactions.hg19, chain.19to38, return.granges = T)
enhancers.to.test <- interactions.hg38[interactions.hg38$med.chicago > 5]
baits.to.test <- extract_promoters_from_interactions(enhancers.to.test)
promoters.to.test <- baits.to.test
dmps.in.enhancer.df <- combine_dmps_with_intervals(dmps.gr, enhancers.to.test)
dmps.in.promoter.df <- combine_dmps_with_intervals(dmps.gr, promoters.to.test)
dmp.counts.in.enhancer.df <- summarize_dmp_counts_in_pchic(dmps.in.enhancer.df, "oe.id", filter.zeros = T)
dmp.counts.in.promoter.df <- summarize_dmp_counts_in_pchic(dmps.in.promoter.df, "bait.id", filter.zeros = T)
genes.w.dm.enhancers <- get_unique_genes_from_baitName_col(dmp.counts.in.enhancer.df)
genes.w.dm.promoters <- get_unique_genes_from_baitName_col(dmp.counts.in.promoter.df)
pchic.summary.stats <- tabulate_pchic_findings(
               enhancers.to.test,
               baits.to.test,
               dmp.counts.in.enhancer.df,
               dmp.counts.in.promoter.df,
               genes.w.dm.enhancers,
               genes.w.dm.promoters
             )
pchic.summary.stats.csv <- my_write_csv(pchic.summary.stats, "/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/MCI_CONTROL/PCHiC-integration-counts.csv")

# compare with DEGs
autosomal.symbols <- get_autosomoal_gene_universe()
diff.exp.data <- clean_differentially_expressed_genes("/media/Data/WGBS/LOAD_MCI/2020-Shigemizu-AD-RNAseq-DEGenes.xlsx",autosomal.symbols)
dmp.enhancer.promoter.integrated <- combine_enhancer_promoter_gwas_diff_exp( dmp.counts.in.enhancer.df,
     dmp.counts.in.promoter.df,
     natgen.symbols,
     diff.exp.data$gene.name
   )
genes.with.dm.enhancer.and.dm.promoter <-  get_common_genes_from_DM_interactions(
                dmps.in.enhancer.df,
                dmps.in.promoter.df)
dm.interactions.summary <- summarize_counts_dm_interactions(dmps.in.enhancer.df, dmps.in.promoter.df)
dm.enhancer.promoter.gene.ont.df <- symbols_to_gene_ontology_routine(genes.with.dm.enhancer.and.dm.promoter)
enhancers.counts <- count_enhancers_with_dmp(dmps.in.enhancer.df)
enhancers.natgen.counts <- count_enhancers_with_dmp(
                subset_interactions_by_diff_exp_data(dmps.in.enhancer.df, diff.exp.data$gene.name))
enhancers.summary <- summarize_interactions_with_dmp(dmps.in.enhancer.df)
dm.enhancers.genes.df <- summarize_dm_interaction_genes(enhancers.summary)

promoters.counts <- count_promoters_with_dmp(dmps.in.promoter.df)
promoters.summary <- summarize_interactions_with_dmp(dmps.in.promoter.df)
dm.promoters.genes.df <- summarize_dm_interaction_genes(promoters.summary)
test.enhancer.enrichment <- test_enrichment_for_dmps(pvals.gr, dmps.gr, enhancers.to.test, B=10000)
test.promoter.enrichment <- test_enrichment_for_dmps(pvals.gr, dmps.gr, promoters.to.test, B=10000)

plot.enhancer.enrichment <- plot_enhancer_enrichment_for_dmps(test.enhancer.enrichment, "/media/Data/test-enhancer-enrichment.pdf")

DE.vs.DME.test <- test_ranks_of_pchic_rnaseq(diff.exp.data, dmps.in.enhancer.df)

#dmps.array.gr <- read_and_cast_madrid_data("./DataReference/madrid_cpgs_list.csv")
#array.gene.symbols <- get_array_genes_with_nearby_dmp(dmps.array.gr)
#wgms.sig.in.array.gr <- subsetByOverlaps(pvals.gr, dmps.array.gr, minoverlap = 2)

#gene.list.by.diff.meth.df <- curate_genes_by_dm_status(dm.genes.df$gene_name, dm.promoters.genes.df$gene.name, dm.enhancers.genes.df$gene.name, diff.exp.data$gene.name, natgen.symbols, array.gene.symbols)

# Export for UCSC
#interactions.for.ucsc <- export_significant_interactions_to_UCSC(enhancers.with.dmp)
#interactions.for.ucsc <- export_significant_interactions_to_UCSC(dmps.in.enhancer.df)
interactions.for.ucsc <- export_significant_interactions_to_UCSC(enhancers.summary)
format_and_write_ucsc_interactions(interactions.for.ucsc,"/media/Data/interactions-with-DMP.hg38.bb")
format_and_write_ucsc_lolly(pvals.gr, "/media/Data/DMPs-lolly.hg38.bb", ALPHA)


# Save and write out results
paper.stats <- get_paper_stats(gene.body.enrichment, DMGENE.ALPHA)
process_and_write_dmps(dmps.gr,
               "Table S1: DMPs.List of Differentially Methylated Positions (DMPs) with coordinates (hg38), effect sizes, and local False-Discovery Rates (lFDRs)",
               "/media/Data/S1-DMPs.xlsx")

process_and_write_nature_genetics_list(
               NG.tally.df,
               natgen.symbols,
               genes.w.dm.enhancers,
               genes.w.dm.promoters,
               "Table S2: List of 75 previously identified genetics risk loci with number of DMPs (if any) within 25,000 nucleotides of gene start/stop",
               "/media/Data/S2-ADRiskLociNumberOfDMPs.xlsx")

process_and_write_DM_genes(
               gene.body.enrichment,
               "Table S3: Differentially methylated genes with coordinates and lFDRs",
               "/media/Data/S3-DMGenes.xlsx")

process_and_write_gene_ontology_terms(
               dm.genes.go.df,
               "Table S4: Gene Ontologies for DM Genes",
               "/media/Data/S4-DMGenes-GeneOntologies.xlsx")

format_and_write_dm_promoters(
               promoters.summary,
               diff.exp.data$gene.name,
               "Table S5: Promoter-enhancer interactions with at least one DMP in promoter",
               "/media/Data/S5-InteractionsWithDMPromoters.xlsx")

format_and_write_dm_enhancers(
               enhancers.summary,
               diff.exp.data$gene.name,
               "Table S6: Promoter-enhancer interactions with at least one DMP in enhancer",
               "/media/Data/S6-InteractionWithDMEnhancers.xlsx")

#my_write_csv(gene.list.by.diff.meth.df, "/media/Data/GenesWithAnyDMFeatures.csv")
