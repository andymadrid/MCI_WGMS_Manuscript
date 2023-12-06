# General purpose functions
# Should be analysis agnostic


# Constants ---------------------------------------------------------------

BLOOD.CELL.TYPES <- c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery",
                      "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8",
                      "nB", "tB")

get_sample_ids_from_dss_inputs <- function(file){
  e <- new.env()
  load(file, envir = e)
  return(rownames(e$design.df))
}


run_chisq_test_from_counts <- function(N.in.a, N.in.b, N.in.both, N.total, labs){
  #
  # Args
  # N.in.a : number of counts in group A
  cell.11 <- N.total - N.in.a - N.in.b + N.in.both
  cell.12 <-  N.in.a - N.in.both
  cell.21 <- N.in.b - N.in.both
  cell.22 <- N.in.both

  X <- matrix(c(cell.11, cell.12, cell.21, cell.22), nrow = 2)
  colnames(X) <- paste0(c("No", "Yes"), labs[1])
  rownames(X) <- paste0(c("No", "Yes"), labs[2])

  chisq.test(X)

  prop.table(X, margin = 1)
}


get_apoe_allele_frequencies <- function(file, sample.ids){
  df <- read_csv(file, show_col_types = F) %>%
    dplyr::filter(sample_id %in% sample.ids) %>%
    dplyr::mutate(APOE.risk.allele = (apoe_e1 == 4) + ( apoe_e2 == 4)) %>%
    dplyr::select(diagnostic_group, sample_id, APOE.risk.allele, apoe_e1, apoe_e2) %>%
    drop_na()
}


test_apoe_allele_frequencies <- function(apoe.df){
  tb <- apoe.df %>%
    group_by(diagnostic_group, APOE.risk.allele) %>%
    summarize(count = n()) %>%
    pivot_wider(names_from = APOE.risk.allele, values_from = count) %>%
    column_to_rownames("diagnostic_group")

  return(list("Table" = tb, "Test" = chisq.test(tb)))
}


download_chain_from_ucsc <- function(dir, file){
  # Assign output name
  ofile.name <- file.path(dir, basename(file))

  # Download
  download.file(file, destfile = ofile.name)

  # Decompress
  R.utils::gunzip(ofile.name, remove = F, overwrite = T)

  # Return a chain object
  rtracklayer::import.chain(stringr::str_remove(ofile.name, "\\.gz"))
}


to_ucsc_format <- function(a, b, c){
  # Return chr1:1111-2222
  if (!str_detect(a, "chr")){a <- paste0("chr", a)}
  paste0(a, ":", b, "-", c)
}

to_ucsc_format_v <- Vectorize(to_ucsc_format)


get_group_ids <- function(master.df, vv){
  sort(unique(master.df$sample_id[master.df$diagnostic_group == vv]))
}


munge_master_df <- function(master.full.df, sample.ids){
  master.full.df %>%
    dplyr::filter(sample_id %in% sample.ids) %>%
    group_by(sample_id) %>%  # one sample shows up multiple times
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::mutate(APOE_risk_allele = (apoe_e1 == 4) + ( apoe_e2 == 4))
}


tabulate_and_test <- function(master.df, var, test.type = "chi"){
  tb <- master.df %>%
    dplyr::count(diagnostic_group, !!rlang::sym(var)) %>%
    pivot_wider(names_from = "diagnostic_group", values_from = "n", values_fill = 0) %>%
    drop_na() %>%
    column_to_rownames(var)

  if (test.type == "chi"){
    test.result <- chisq.test(tb)
  } else if (test.type == "fisher"){
    test.result <- fisher.test(tb)
  }

  return(list("Table" = tb,
              "Test" = test.result))
}


stringify_mean_sd <- function(xx){
  paste0(round(mean(xx), 2), " (", round(sd(xx), 2), ")")
}


test_continuous <- function(master.df, var){

  if (!(var %in% colnames(master.df))){
    warning("Var not in master.df")
  }

  xx.control <- master.df %>%
    dplyr::filter(diagnostic_group == "CONTROL") %>%
    dplyr::pull(var) %>% as.numeric()

  xx.mci <- master.df %>%
    dplyr::filter(diagnostic_group == "MCI") %>%
    dplyr::pull(var) %>% as.numeric()

  list("Control: " = stringify_mean_sd(xx.control),
       "MCI: " = stringify_mean_sd(xx.mci),
       "Test: " = t.test(xx.control, xx.mci))
}

test_continuous_three_group <- function(master.df, var){

  if (!(var %in% colnames(master.df))){
    warning("Var not in master.df")
  }
  if (var == "education") {
    master.df <- master.df %>%
    dplyr::filter(!is.na(education))
  }
  xx.control <- master.df %>%
    dplyr::filter(diagnostic_group == "CONTROL") %>%
    dplyr::pull(var) %>% as.numeric()

  xx.load <- master.df %>%
    dplyr::filter(diagnostic_group == "LOAD") %>%
    dplyr::pull(var) %>% as.numeric()

  xx.mci <- master.df %>%
    dplyr::filter(diagnostic_group == "MCI") %>%
    dplyr::pull(var) %>% as.numeric()
if (var == "bmi") {
  list("Control: " = stringify_mean_sd(xx.control),
       "MCI: " = stringify_mean_sd(xx.mci),
       "LOAD: " = stringify_mean_sd(xx.load),
       "Test: " = pairwise.t.test(master.df$bmi, master.df$diagnostic_group, p.adjust.method="none"))
}
else if (var == "age_at_visit") {
  list("Control: " = stringify_mean_sd(xx.control),
       "MCI: " = stringify_mean_sd(xx.mci),
       "LOAD: " = stringify_mean_sd(xx.load),
       "Test: " = pairwise.t.test(master.df$age_at_visit, master.df$diagnostic_group, p.adjust.method="none"))
}
else if (var == "education") {
  list("Control: " = stringify_mean_sd(xx.control),
       "MCI: " = stringify_mean_sd(xx.mci),
       "LOAD: " = stringify_mean_sd(xx.load),
       "Test: " = pairwise.t.test(master.df$education, master.df$diagnostic_group, p.adjust.method="none"))
}
}
# Missingness quanitifcation ----------------------------------------------


read_and_tally_subroutine <- function(files){
  all <- do.call(rbind,
                 lapply(X = files,
                        FUN = fread,
                        select = c("chromStart", "coverage", "sample")))

  wide.df <- pivot_wider(all,
                         "chromStart",
                         names_from = "sample",
                         values_from = "coverage")

  out <- data.frame(chr = "chr1",
             start = wide.df$chromStart,
             end = wide.df$chromStart + 2,
             N.missing = rowSums(is.na(wide.df)))

  print("Read and pivoted one chunk...")
  return(out)
}

read_and_join_m_cov_routine <- function(dir, sample.ids, chunk.size){
  all.files <- file.path(dir, paste0("chr1.", sample.ids, ".bed"))

  all.files.split <- split(all.files, ceiling(seq_along(all.files) / chunk.size))

  counts.df <- do.call(rbind, lapply(X = all.files.split, read_and_tally_subroutine)) %>%
    group_by(chr, start, end) %>%
    summarize(PercentMissing = sum(N.missing) / length(sample.ids))

  return(counts.df)
}


# Statistical functions ---------------------------------------------------

compute_empirical_p <- function(xx, test.value){
  positive.test.value <- abs(test.value)

  p1 <- (xx >= test.value) / length(xx)
  p2 <- (xx <= test.value) / length(xx)

  p <- 2 * min(p1, p2)
  return(p)
}










subset_dmps <- function(data, lfdr.cut=0.05){
  data %>%
    dplyr::filter(lfdr <= lfdr.cut) %>%
#    dplyr::filter (abs(pi.diff.mci.ctrl) >= DMALPHA)
    dplyr::filter (abs(pi.diff) >= DMALPHA)
}


to_granges <- function(data){
  out <- data %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns=T,
      starts.in.df.are.0based = T
    )

  seqlevelsStyle(out) <- "UCSC"
  out
}





# Get public data ---------------------------------------------------------

clean_differentially_expressed_genes <- function(file, autosomal.symbols){
  readxl::read_xlsx(file, skip=1) %>%
    dplyr::rename(gene.name = "gene name",
                  gene.type = "gene type",
                  entrez.id = "Entrez Gene ID",
                  gene.id = "Ensemble ID",
                  expression.pval = "P-value") %>%
    dplyr::filter(gene.name %in% autosomal.symbols)
}

# Essentially constants ---------------------------------------------------


get_hyper_hypo_colors <- function(rgb = F){
  # Hyper and hypo colors, may be used outside volcano
  # so we'll keep it in the R/functions.R file

  if (rgb){
    return(list(hyper = "239,192,0", hypo = "0,119,154"))
  } else {
    return(list(hyper = "#0073C2FF", hypo = "#EFC000FF"))
  }
}



# Manipulate GRanges ------------------------------------------------------

get_ensdb <- function(){
  # Make it super standard
  return(EnsDb.Hsapiens.v86)
}

get_autosomoal_gene_universe <- function(protein_coding=F){
  # if protein_coding = T, return only protein_coding genes

  edb <- get_ensdb()

  all.genes <- genes(edb)
  autosomes <- all.genes[seqnames(all.genes) %in% 1:22]

  if (protein_coding){
    autosomes$symbol[autosomes$gene_biotype == "protein_coding"]
  } else {
    autosomes$symbol
  }
}

expand_genes <- function(gr, upstream, downstream) {
  # https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  strand_is_minus = strand(gr) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(gr)[on_plus] = start(gr)[on_plus] - upstream
  start(gr)[on_minus] = start(gr)[on_minus] - downstream
  end(gr)[on_plus] = end(gr)[on_plus] + downstream
  end(gr)[on_minus] = end(gr)[on_minus] + upstream
  return(gr)
}

get_gene_bodies <- function(upstream, downstream, autosomes_only = T, protein_coding_only = T){
  db <- get_ensdb()
  genes <- genes(db)

  # Respective indices
  protein.coding.ix <- genes$gene_biotype == "protein_coding"
  autosome.ix <- as.vector(seqnames(genes) %in% 1:22)

  if (autosomes_only & protein_coding_only){
    out <- expand_genes(genes[protein.coding.ix & autosome.ix], upstream, downstream)
  } else if (autosomes_only & !protein_coding_only){
    out <- expand_genes(genes[autosome.ix], upstream, downstream)
  } else if (!autosomes_only & protein_coding_only){
    out <- expand_genes(genes[protein.coding.ix], upstream, downstream)
  } else {
    out <- expand_genes(genes, upstream, downstream )
  }

  seqlevelsStyle(out) <- "NCBI"
#   seqlevelsStyle(out) <- "UCSC"
 out
}

ids_to_symbols <- function(vv){
  genes.ref <- ensembldb::genes(get_ensdb())
  genes.ref[match(vv, genes.ref$gene_id )]$symbol
}

symbols_to_ids <- function(vv){
  genes.ref <- ensembldb::genes(get_ensdb())
  ix <- match(vv, genes.ref$symbol)
  genes.ref[ix[!is.na(ix)]]$gene_id
}

get_protein_coding_promoters <- function(upstream, downstream){
  db <- get_ensdb()
  promoters <- promoters(db, upstream, downstream)
  promoters <- promoters[promoters$tx_biotype == "protein_coding"]
  promoters$gene_name <- ids_to_symbols(promoters$gene_id)
  promoters
}


make_df_from_two_overlapping_granges <- function(range.1, range.2){
  # Common seqlevels
  seqlevelsStyle(range.1) <- "NCBI"
  seqlevelsStyle(range.2) <- "NCBI"

  # Overlap and get indices
  overlaps <- findOverlaps(range.1, range.2)
  ix.1 <- queryHits(overlaps)
  ix.2 <- subjectHits(overlaps)

  # We'll keep the seqnames for range.1, not range.2 (otherwise duplicate column names)
  left <- data.frame(range.1[ix.1, ]) %>%
    dplyr::select(-c(width, strand)) %>%
    dplyr::rename(chr = seqnames)

  right <- data.frame(range.2)[ix.2, -1:-3] %>%
    dplyr::select(-width)

  cbind(left, right)
}

combine_missingness_with_pvals <- function(missing.chr1.load.df,
                                           missing.chr1.control.df,
                                           pvals.gr){

  colnames(missing.chr1.load.df)[4] <- "AD"
  colnames(missing.chr1.control.df)[4] <- "NoAD"

  miss.df <- full_join(missing.chr1.load.df,
                               missing.chr1.control.df,
                               by = c("chr", "start", "end"))

  miss.gr <- makeGRangesFromDataFrame(miss.df, keep.extra.columns = T)
  make_df_from_two_overlapping_granges(miss.gr, pvals.gr)
}

plot_missingness_hexbin <- function(miss.data, file){
  z <- miss.data %>%
    pivot_longer(contains("AD"), values_to = "Percent missing", names_to = "Group") %>%
    ggplot(aes(x = `Percent missing`, y = -log10(lfdr))) +
      facet_wrap(.~Group) +
      stat_binhex() +
      theme_minimal() +
      scale_fill_viridis(option = "B", trans = "log", breaks = 10^(1:5)) +
      ylab(expression(-log[10](lFDR))) +
      geom_hline(yintercept = -log10(0.05), color = "grey40") +
    ggtitle("Percent missing values before imputation") +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  cowplot::save_plot(file, z)
  return(file)
}


# Harmonic P-value routine ------------------------------------------------

harmonic_pvalue_routine <- function(loci.gr, features.gr, alpha){
  df <- make_df_from_two_overlapping_granges(loci.gr, features.gr)

  hmp.df <- df %>%
    group_by(gene_name) %>%
    dplyr::summarize(
      N.CpGs = n(),
#      N.DMPs = sum(lfdr < alpha & abs(pi.diff.mci.ctrl) > DMALPHA),
       N.DMPs = sum(lfdr < alpha & abs(pi.diff) > DMALPHA),
     HarmonicMeanPval = harmonicmeanp::p.hmp(pval, L = N.CpGs)
    ) %>%
    dplyr::mutate(HarmonicMeanPval = pmin(HarmonicMeanPval, 1)) %>%
    dplyr::mutate(HarmonicMeanPval = pmax(HarmonicMeanPval, 0))

  # FDR calculation
  fdr.out <- fdrtool::fdrtool(x = hmp.df$HarmonicMeanPval, statistic = "pvalue", plot = F)

  # Pack together
  hmp.df$lfdr <- fdr.out$lfdr
  hmp.df$qval <- fdr.out$qval

  # Output
  hmp.df
}

harmonic_pvalue_routine_hyper_hypo <- function(loci.gr, features.gr, alpha){
  df <- make_df_from_two_overlapping_granges(loci.gr, features.gr)

  # Hyper run
  hmp.df <- df %>%
    group_by(gene_name) %>%
    dplyr::summarize(
      N.CpGs = n(),
      N.DMPs = sum(lfdr < alpha & pi.diff > DMALPHA),
      HarmonicMeanPval = harmonicmeanp::p.hmp(pval, L = N.CpGs)
    ) %>%
    dplyr::mutate(HarmonicMeanPval = pmin(HarmonicMeanPval, 1)) %>%
    dplyr::mutate(HarmonicMeanPval = pmax(HarmonicMeanPval, 0))

  # FDR calculation
  fdr.out <- fdrtool::fdrtool(x = hmp.df$HarmonicMeanPval, statistic = "pvalue", plot = F)

  # Pack together
  hmp.df$lfdr <- fdr.out$lfdr
  hmp.df$qval <- fdr.out$qval

  # Output
#  hmp.df
  hmp.df.hyper <- hmp.df

  # Hypo run
  hmp.df <- df %>%
    group_by(gene_name) %>%
    dplyr::summarize(
      N.CpGs = n(),
      N.DMPs = sum(lfdr < alpha & pi.diff < (DMALPHA*-1)),
      HarmonicMeanPval = harmonicmeanp::p.hmp(pval, L = N.CpGs)
    ) %>%
    dplyr::mutate(HarmonicMeanPval = pmin(HarmonicMeanPval, 1)) %>%
    dplyr::mutate(HarmonicMeanPval = pmax(HarmonicMeanPval, 0))

  # FDR calculation
  fdr.out <- fdrtool::fdrtool(x = hmp.df$HarmonicMeanPval, statistic = "pvalue", plot = F)

  # Pack together
  hmp.df$lfdr <- fdr.out$lfdr
  hmp.df$qval <- fdr.out$qval

  # Output
#  hmp.df
  hmp.df.hypo <- hmp.df
  hmp.df <- list("hyper" = hmp.df.hyper, "hypo" = hmp.df.hypo)
  hmp.df
}


# -------------------------------------------------------------------------

filter_by_natgen<- function(data, genes){
  data %>%
    dplyr::filter(gene_name %in% genes) %>%
    arrange(lfdr)
}



# Gene ontology pipeline --------------------------------------------------

run_gene_ontology <- function(genes){
  clusterProfiler::enrichGO(
    genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "fdr"
  )
}

go_output_to_df <- function(go.out){
  # Constant to define mapping from abbreviated terms to full ones used in plots
  MAPPING <- data.frame(
    ONTOLOGY = c("MF", "CC", "BP"),
    Ontology = c("Molecular\nfunction", "Cellular\ncomponent", "Biological\nprocess")
  )

  go.df <- data.frame(go.out)

  # Get the names right
  dplyr::left_join(go.df, MAPPING, by = "ONTOLOGY")
}


plot_go_barchart <- function(go.df, n=25, ofile){
  height = 7
  width = 14
  # Get top 25 by p-value, then arrang by gene set size
  subdata <- head(dplyr::arrange(go.df, p.adjust), n) %>%
    dplyr::arrange(-p.adjust) %>%
    dplyr::mutate(Count = as.numeric(Count))

  subdata$Description <- factor(subdata$Description, levels = subdata$Description)

  pdf(ofile,height=7,width=14)
  p <- ggplot(data = subdata,
         aes(x = Description,
#             y = Count,
             y = -log10(p.adjust),
             fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_meAD() +
    xlab("") +
#    ylab("Number of genes") +
    ylab("-log10(FDR)") +
    labs(legend = "") +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "white", color = "white")) +
    scale_fill_manual(values = c("firebrick2", "mediumblue", "darkgoldenrod2"))
  print(p)
  dev.off()
}




symbols_df_to_go_df_routine <- function(file){

  # Assuming we read from file
  df <- read_csv(file, show_col_types = F)
  gene.ids <- symbols_to_ids(df$gene_name)

  # Gene ontology and data munging
  go.out <- run_gene_ontology(gene.ids)
  go_output_to_df(go.out)

}


convert_gene_ontology_ids_to_symbols <- function(DMGenes.go.df){
  # Return the same data frame but with gene symbols instead of gene ontologies
  DMGenes.go.df %>%
    mutate(gene.ids = str_split(geneID, "/")) %>%
    unnest(gene.ids) %>%
    mutate(gene.symbols = ids_to_symbols(gene.ids)) %>%
    group_by(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust) %>%
    summarize(GeneSymbols = paste(gene.symbols, collapse = ";"))
}


symbols_to_gene_ontology_routine <- function(symbols){
  symbols_to_ids(symbols) %>%
    run_gene_ontology %>%
    go_output_to_df %>%
    convert_gene_ontology_ids_to_symbols %>%
    tidyr::separate(GeneRatio, c("Count", "BGCount"), remove = F) %>%
    dplyr::rename(Ontology = "ONTOLOGY")
}


# Export UCSC Genome Browser ----------------------------------------------

get_ucsc_seqlengths <- function(){
  good.chr <- paste0("chr", 1:22)

  db <- get_ensdb()
#  seqlevelsStyle(db) <- "NCBI"
  seqlevelsStyle(db) <- "UCSC"

  # Usually pings to say "[some weird scaffold] didn't map"
  suppressWarnings(seqlengths(db)[good.chr])
}


#TODO: refactor
format_and_write_ucsc_lolly <- function(data.gr, file, lfdr.cut, keep.nth=25){
  # data.gr is pvals.gr or something similar
  # alpha.cut is the level of significance
  # keep.nth is how much to thin by (keep every 25th non-significant point)

  colors <- get_hyper_hypo_colors(rgb = T)

  # Thin the points
  keepix <- c(
    which(data.gr$lfdr <= lfdr.cut),
    which(data.gr$lfdr > lfdr.cut)[c( rep(FALSE, keep.nth-1), TRUE)]
  )


  # Thin data
  out <- data.gr[keepix, ]

  # Manipulate the metadata columns to get in UCSC lolly format
  mdata.cleaned <- mcols(out) %>%
    as.data.frame() %>%
    dplyr::transmute(name = ".",
                     score = round(y),
                     thickStart = start(out),
                     thickEnd = end(out),
#                     color = ifelse(pi.diff < 0, colors$hypo, colors$hyper),
#                     color = ifelse(pi.diff.mci.ctrl < 0, colors$hypo, colors$hyper),
                      color = ifelse(pi.diff < 0, colors$hypo, colors$hyper),
                    lollySize = ifelse(lfdr < lfdr.cut, 4, 1))

  mdata.cleaned$color[mdata.cleaned$lollySize == 1] <- "220,220,220"

  mcols(out) <- mdata.cleaned
  genome(out) <- "hg38"
  seqlengths(out) <- get_ucsc_seqlengths()

  # Write and return
  rtracklayer::export.bb(out, file)
  return(file)
}


liftover_wrapper <- function(gr, chain){
  seqlevelsStyle(gr) <- "UCSC"

  unlist(rtracklayer::liftOver(gr, chain)) %>%
    as.data.frame()

}

filter_by_gene_symbols <- function(genes, inclusion.genes){
  out <- genes[genes$symbol %in% inclusion.genes]
  out
}


tally_dmps_in_genes <- function(genes.with.dmp){
  as.data.frame(genes.with.dmp) %>%
    group_by(symbol) %>%
    summarize(N.DMPs = n())
}



get_array_genes_with_nearby_dmp <- function(dmps.array.gr){
  sort(unique(unlist(str_split(dmps.array.gr$gene.symbols, ";"))))
}

compare_with_madrid_paper <- function(dmps.array.gr, dm.genes.df){
  array.genes <- get_array_genes_with_nearby_dmp(dmps.array.gr)
  wgms.genes <- dm.genes.df$gene_name

  zz <- intersect(array.genes, wgms.genes)

  return(zz)
}


create_wb_with_description <- function(desc){
  cc <- 1:6

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, 1)
  openxlsx::writeData(wb, sheet = 1,
                      startCol = 1,
                      startRow = 1,
                      desc,
                      headerStyle = )


  openxlsx::addStyle(wb, sheet = 1, cols = cc, rows = 1, createStyle(textDecoration = "Bold"))

  return(wb)
}

append_data_to_workbook <- function(wb, data){



  setColWidths(wb, sheet = 1, cols = 1:ncol(data), widths = pmax(nchar(colnames(data)) * 1.2, 8) )
  mergeCells(wb, sheet=1, cols = 1:ncol(data), rows=1)


  openxlsx::writeData(wb = wb, sheet=1, x = data, startRow = 2,
                      headerStyle = createStyle(textDecoration = "Bold",
                                                halign = "center"))
  return(wb)
}

tally_and_write_DMGenes_with_DMP_counts <- function(
    gene.body.enrichment,
    lfdr.cut,
    desc,
    file){

  data <- gene.body.enrichment %>%
    dplyr::filter(lfdr < lfdr.cut)

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}


run_length_embed_dmps <- function(dmps.data, chr, dist){
  tmp <- dmps.data %>%
    dplyr::filter(chr == chr) %>%
    arrange(start)

  rle(diff(tmp$start) < dist)
}


get_N_not_in_set <- function(dmps.gr, features.gr){
  seqlevelsStyle(dmps.gr) <- "NCBI"
  seqlevelsStyle(features.gr) <- "NCBI"

  out <- subsetByOverlaps(x = dmps.gr, ranges = features.gr, invert = TRUE)
  return(length(out))
}

get_N_in_set <- function(dmps.gr, features.gr){
  seqlevelsStyle(dmps.gr) <- "NCBI"
  seqlevelsStyle(features.gr) <- "NCBI"

  out <- subsetByOverlaps(x = dmps.gr, ranges = features.gr, invert = FALSE)
  return(length(out))
}


tally_dmps_in_out_genes <- function(dmps.gr, features.gr){
  return(list(
    "N DMPs  in set:" = get_N_in_set(dmps.gr, features.gr),
    "N DMPs NOT in set:" = get_N_not_in_set(dmps.gr, features.gr)
  ))
}

#END



ucsc.order.cols <- c("chrom", "chromStart", "chromEnd",
                     "name", "score", "value", "exp", "color",
                     "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand",
                     "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand")

# unnest_interactions <- function(interactions.gr){
#   interactions.gr %>%
#     dplyr::mutate(gene.name = str_split(baitName, ";")) %>%
#     tidyr::unnest(gene.name)
# }


clean_interactions_data <- function(file){
  # Read as text file
  promoter_capture.data <- fread(file)

  # Filter out the following:
  # 1) oeName == "." indicates a bait-to-bait interaction
  # 2) baitName == NA indicates unknown promoter
  # 3) dist == NA indicates trans-chromosomal interaction
  # 4) Sex chromosomes (don't include for now)
  promoter_capture.cleaned <- promoter_capture.data %>%
    dplyr::filter(oeName == ".",
                  !is.na(baitName),
                  # !is.na(dist),
                  !(baitChr %in% c("X", "Y"))) %>%
    # Add identifiers for each bait and other end
    dplyr::mutate(bait.id = to_ucsc_format_v(baitChr, baitStart, baitEnd),
                  oe.id = to_ucsc_format_v(oeChr, oeStart, oeEnd)) %>%
    # Drop unneeded columns
    dplyr::select(-c(clusterID, clusterPostProb, oeName, oeID, baitID))

  # Get some summary statistics
  promoter_capture.cleaned$min.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=min)
  promoter_capture.cleaned$med.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=median)
  promoter_capture.cleaned$mean.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=mean)
  promoter_capture.cleaned$max.chicago <- apply(X=promoter_capture.cleaned[, ..BLOOD.CELL.TYPES], MARGIN=1, FUN=max)

  # Return this big fella
  promoter_capture.cleaned
}



make_granges_with_common_field_prefix <- function(gr, prefix=""){
  # Helper for interactions data. Basically allows you to specify the predix
  # so you can say "bait" and this function will figure out that
  # the output GRanges nees baitChr, baitStart, baitEnd

  makeGRangesFromDataFrame(gr,
                           seqnames.field =  paste0(prefix, "Chr"),
                           start.field = paste0(prefix, "Start"),
                           end.field = paste0(prefix, "End"),
                           keep.extra.columns = T)
}



drop_granges_columns <- function(data){
  # Remove these columns
  dplyr::select(data, -c(seqnames, start, end, width, strand))
}


lift_promoter_capture_data_to_hg38 <- function(interactions.hg19, chain, return.granges=F){
  # Need an id to re-combine later
  interactions.hg19$interaction.id <- 1:nrow(interactions.hg19)

  #TODO: how to reconcile differences here??
  baits.hg19 <- make_granges_with_common_field_prefix(interactions.hg19, prefix="bait")

  # Convert to UCSC style and liftOver
  baits.hg38 <-liftover_wrapper(baits.hg19, chain)
  baits.hg38 <- baits.hg38[which(!duplicated(baits.hg38$interaction.id)==TRUE),]
  baits.hg38$bait.id <- paste0(baits.hg38$seqnames,":",baits.hg38$start,"-",baits.hg38$end)
  baits.hg38 <- baits.hg38 %>%
	drop_granges_columns() %>%
    distinct()
  baits.hg38 <- baits.hg38[,-which(colnames(baits.hg38)=="oe.id")]
  # Other ends
  other_ends.hg19 <- make_granges_with_common_field_prefix(interactions.hg19, prefix="oe")

  # Convert to UCSC style and liftOver
  other_ends.hg38 <- liftover_wrapper(other_ends.hg19, chain)
  other_ends.hg38 <- other_ends.hg38[which(!duplicated(other_ends.hg38$interaction.id)==TRUE),]
  other_ends.hg38$oe.id <- paste0(other_ends.hg38$seqnames,":",other_ends.hg38$start,"-",other_ends.hg38$end)
  other_ends.hg38$baitChr <- other_ends.hg38$seqnames
  other_ends.hg38$baitStart <- other_ends.hg38$start
  other_ends.hg38$baitEnd <- other_ends.hg38$end
  other_ends.hg38 <- other_ends.hg38 %>%
    drop_granges_columns() %>%
    dplyr::select(c("interaction.id", "baitChr", "baitStart", "baitEnd", "oe.id")) %>%
    distinct()

  # Merge
  output <- dplyr::inner_join(baits.hg38, other_ends.hg38, by = "interaction.id")
  output$oeChr <- output$baitChr
  output$oeStart <- output$baitStart
  output$oeEnd <- output$baitEnd
  output$oeChr <- gsub("chr","",output$oeChr)
  output$baitChr <- gsub("chr","",output$baitChr)
	
  # Change bait information
  bait.info <- as.data.frame(output$bait.id)
  bait.info <- as.data.frame(str_split_fixed(bait.info[,1],":",2))
  bait.info.coords <- as.data.frame(str_split_fixed(bait.info[,2],"-",2))
  bait.info <- cbind(bait.info[,1],bait.info.coords)
  colnames(bait.info) <- c("chr","start","end")	
  output$baitChr <- bait.info$chr
  output$baitStart <- bait.info$start
  output$baitEnd <- bait.info$end
  output$baitChr <- gsub("chr","",output$baitChr)

  if (nrow(output) == length(unique(output$interaction.id))){
    if (return.granges){
      return(make_granges_with_common_field_prefix(output, prefix = "oe"))
    } else {
      return(output)
    }
  } else {
    warning("N row of interactions is not the same as unique interactions")
  }
}

#lift_promoter_capture_data_to_hg38 <- function(interactions.hg19, chain, return.granges=F){
  # Need an id to re-combine later
#  interactions.hg19$interaction.id <- 1:nrow(interactions.hg19)

  #TODO: how to reconcile differences here??
#  baits.hg19 <- make_granges_with_common_field_prefix(interactions.hg19, prefix="bait")

  # Convert to UCSC style and liftOver
#  baits.hg38 <-liftover_wrapper(baits.hg19, chain) %>%
#    drop_granges_columns() %>%
#    distinct()

  # Other ends
#  other_ends.hg19 <- make_granges_with_common_field_prefix(interactions.hg19, prefix="oe")

  # Convert to UCSC style and liftOver
#  other_ends.hg38 <- liftover_wrapper(other_ends.hg19, chain) %>%
#    drop_granges_columns() %>%
#    dplyr::select(c("interaction.id", "baitChr", "baitStart", "baitEnd")) %>%
#    distinct()

  # Merge
#  output <- dplyr::inner_join(baits.hg38, other_ends.hg38, by = "interaction.id")

#  if (nrow(output) == length(unique(output$interaction.id))){
#    if (return.granges){
#      return(make_granges_with_common_field_prefix(output, prefix = "oe"))
#    } else {
#      return(output)
#    }
#  } else {
#    warning("N row of interactions is not the same as unique interactions")
#  }
#}



subset_overlaps_return_unique <- function(features.gr, dmps.gr){
  seqlevelsStyle(dmps.gr) <- "UCSC"
  seqlevelsStyle(features.gr) <- "UCSC"

  subsetByOverlaps(features.gr, dmps.gr) %>%
    unique()
}




count_unique_overlaps <- function(dmps.gr, features.gr){
  # Helper function for null distribution (enrichment analysis)
  seqlevelsStyle(dmps.gr) <- "UCSC"
  seqlevelsStyle(features.gr) <- "UCSC"

  overlaps <- findOverlaps(dmps.gr, features.gr)

  # Tally
  return(list(
    "N.features" = length(unique(subjectHits(overlaps))),
    "N.dmps" = length(unique(queryHits(overlaps)))
  ))
}



# Enrichment Analysis -----------------------------------------------------

test_enrichment_for_dmps <- function(null.space.gr, dmps.gr, features.gr, B=10000){

  N.simulated.dmps <- length(dmps.gr)

  # All possible DMPs in simulation
  null.ix <- 1:length(null.space.gr)

  N.enhancers <- rep(NA, B)
  N.dmps <- rep(NA, B)

  # Run the simulation
  for (b in 1:B){
    sim.ix <- sample(null.ix, size = N.simulated.dmps, replace = F)
    my.counts <- count_unique_overlaps(null.space.gr[sim.ix, ], features.gr)

    # Tally
    N.enhancers[b] <- my.counts$N.features
    N.dmps[b] <- my.counts$N.dmps
  }
  return(list("simulated.df" = data.frame(N.enhancers, N.dmps),
              "test.stat" = count_unique_overlaps(dmps.gr, features.gr)))
}

plot_enhancer_enrichment_for_dmps <- function(enrichment.result, file) {
  xx <- enrichment.result$simulated.df$N.dmps
  test.value <- enrichment.result$test.stat$N.dmps
  pval <- compute_empirical_p(xx, test.value)

  # Clean up
  if (pval == 0){
    pval.str <- "< 0.001"
  } else {
    pval.str <- round(pval, digits = 2)
  }

  z <- enrichment.result$simulated.df %>%
    ggplot(aes(x = N.dmps, y = after_stat(density))) +
    geom_histogram(bins = 25) +
    theme_minimal() +
    geom_vline(xintercept = enrichment.result$test.stat$N.dmps, color = "red") +
    xlab("Number of unique DMPs residing in feature") +
    ylab("Density") +
    labs(caption = paste0("Two-sided p-value: ", pval.str)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

#  cowplot::save_plot(filename = file, z)
  pdf(file)
    print(z)
    dev.off()
}



combine_dmps_with_intervals <- function(dmps.gr, interactions.gr){
  # Overlap DMPs and interactions, and return a data frame.
  # There will be one row for each dmp by interaction (so some duplication if an enhancer
  # is overlapped by multiple DMPs)
  make_df_from_two_overlapping_granges(dmps.gr, interactions.gr) %>%
    group_by(interaction.id) %>%
    dplyr::mutate(
#      median.pi.diff = median(pi.diff.mci.ctrl),
       median.pi.diff = median(pi.diff),
     k.dmps = n()) %>%
    ungroup()
}


get_unique_genes_from_baitName_col <- function(interactions.data){
  sort(unique(unlist(str_split(interactions.data$baitName, ";"))))
}

summarize_interactions_with_dmp <- function(interactions.with.dmps) {
  interactions.with.dmps %>%
    # Drop DMP related stuff
    dplyr::select(-c(chr, start, end, stat, pval, pval.Wald, y, strand,
  #                   diagnostic_group_coded, pi.diff, lfdr)) %>%
#                     diagnostic_groupLOAD, pi.diff.mci.ctrl, lfdr)) %>%
                      diagnostic_groupLOAD, pi.diff, lfdr)) %>%
 group_by(interaction.id) %>%
    dplyr::slice(1)
}


subset_interactions_by_diff_exp_data <- function(interactions, genes.of.interest){
  interactions %>%
    dplyr::mutate(gene.name = str_split(baitName, ";")) %>%
    tidyr::unnest(gene.name) %>%
    dplyr::filter(gene.name %in% genes.of.interest) %>%
    dplyr::select(-gene.name) %>%
    distinct()
}


count_enhancers_with_dmp <- function(enhancers.with.dmp){
  cpg.ids <- paste0(enhancers.with.dmp$chr, ":", enhancers.with.dmp$start)
  enhancer.ids <- enhancers.with.dmp$oe.id
  interaction.ids <- enhancers.with.dmp$interaction.id

  genes <- enhancers.with.dmp %>%
    unnest_interactions_by_gene() %>%
    dplyr::pull(gene.name)

  return(list(
    paste0("N unique DMPs: ", length(unique(cpg.ids))),
    paste0("N unique enhancers: ", length(unique(enhancer.ids))),
    paste0("N interactions: ", length(unique(interaction.ids))),
    paste0("N unique gene: ", length(unique(genes)))
  ))
}


count_promoters_with_dmp <- function(promoters.with.dmp){
  cpg.ids <- paste0(promoters.with.dmp$chr, ":", promoters.with.dmp$start)
  promoter.ids <- promoters.with.dmp$bait.id
  interaction.ids <- promoters.with.dmp$interaction.id

  genes <- promoters.with.dmp %>%
    unnest_interactions_by_gene() %>%
    dplyr::pull(gene.name)

  return(list(
    paste0("N unique DMPs: ", length(unique(cpg.ids))),
    paste0("N unique promoters: ", length(unique(promoter.ids))),
    paste0("N interactions: ", length(unique(interaction.ids))),
    paste0("N unique gene: ", length(unique(genes)))
  ))
}



summarize_counts_dm_interactions <- function(enhancers.with.dmp, promoters.with.dmp){

  symbols <- get_common_genes_from_DM_interactions(enhancers.with.dmp, promoters.with.dmp)

  N.enhancer.genes <- enhancers.with.dmp %>%
    unnest_interactions_by_gene() %>%
    dplyr::pull(gene.name) %>%
    unique() %>% length()

  N.promoter.genes <- promoters.with.dmp %>%
    unnest_interactions_by_gene() %>%
    dplyr::pull(gene.name) %>%
    unique() %>% length()

  return(
    list(
      "N interactions w DM enhancer: " = length(unique(enhancers.with.dmp$interaction.id)),
      "N interactions w DM promoter: " = length(unique(promoters.with.dmp$interaction.id)),
      "N interactions w both: " = length(intersect(enhancers.with.dmp$interaction.id,
                                                   promoters.with.dmp$interaction.id)),
      "N genes w DM enhancer" = N.enhancer.genes,
      "N genes w DM promoter" = N.promoter.genes,
      "N genes w both:" = length(symbols)
      )
    )
}


intersect_semicolon_wrapper <- function(xx, yy){

  xx.split <- unlist(str_split(xx, ";"))

  out <- intersect(xx.split, yy)

  if (length(out) == 0){
    return("None")
  } else {
    paste0(out, collapse = ";")
  }
}


combine_enhancer_promoter_gwas_diff_exp <- function(
    dmp.counts.in.enhancer.df,
    dmp.counts.in.promoter.df,
    gwas.genes,
    diff.exp.genes){

  full_join(dmp.counts.in.enhancer.df,
            dmp.counts.in.promoter.df,
            by = c("interaction.id", "bait.id", "oe.id", "baitName")) %>%
    replace_na(list(n.dmps.x = 0, n.dmps.y = 0)) %>%
    dplyr::mutate(RiskGenes = intersect_semicolon_wrapper(baitName, gwas.genes),
                  DiffExpGenes = intersect_semicolon_wrapper(baitName, diff.exp.genes)) %>%
    dplyr::rename(EnhancerLocus = oe.id,
                  PromoterLocus = bait.id,
                  BaitName = baitName,
                  NumDMPsInEnhancer = n.dmps.x,
                  NumDMPsInPromoter = n.dmps.y,
                  )
}


summarize_dmp_counts_in_pchic <- function(dmps.in.promoter.df, feature.name, filter.zeros){

if (!(feature.name %in% c("oe.id", "bait.id"))){
  stop("Invalid fature.name in summarize_dmp_counts")
}

  out <- dmps.in.promoter.df %>%
    dplyr::mutate(cpg.id = paste0(chr, ":", start)) %>%
    group_by(interaction.id,
             oe.id,
             bait.id,
             baitName) %>%
    summarize(n.dmps = n_distinct(cpg.id))

  if (!filter.zeros){
    out
  } else {
    dplyr::filter(out, n.dmps > 0)
  }
}

tabulate_pchic_findings <- function(
    enhancers.to.test,
    baits.to.test,
    dmp.counts.in.enhancer.df,
    dmp.counts.in.promoter.df,
    genes.w.dm.enhancers,
    genes.w.dm.promoters
    ){

  N.row <- 10
  df <- data.frame(description = rep(NA, N.row), count = rep(NA, N.row))

  # Number of interactions and whatnot tested
  df[1, ] <- c("Interactions included in analysis", n_distinct(enhancers.to.test$interaction.id))
  df[2, ] <- c("Enhancers included in analysis", n_distinct(enhancers.to.test$oe.id))
  df[3, ] <- c("Promoters included in analysis", n_distinct(baits.to.test$bait.id))


  # First, the counts of enhancers
  df[4, ] <- c("Interactions w/ DM enhancer (1+ DMP)",
             n_distinct(dmp.counts.in.enhancer.df$interaction.id))
  df[5, ] <- c("Unique enhancers w/ 1+ DMP",
               n_distinct(dmp.counts.in.enhancer.df$oe.id))

  # Then the promoters
  df[6, ] <- c("Interactions w/ DM promoter (1+ DMP)",
               n_distinct(dmp.counts.in.promoter.df$interaction.id))
  df[7, ] <- c("Unique promoters w/ 1+ DMP",
               n_distinct(dmp.counts.in.promoter.df$bait.id))

  # Now the genes
  df[8, ] <- c("Genes with 1+ DMP in 1+ enhancer", n_distinct(genes.w.dm.enhancers))
  df[9, ] <- c("Genes with 1+ DMP in 1+ promoter", n_distinct(genes.w.dm.promoters))

  df[10, ] <- c("Genes w/ 1+ DMP in both enhancer and promoter",
                length(intersect(genes.w.dm.enhancers, genes.w.dm.promoters)))

  return(df)

}



summarize_dm_interaction_genes <- function(interactions.summary){
  interactions.summary %>%
    dplyr::mutate(gene.name = str_split(baitName, ";")) %>%
    unnest(gene.name) %>%
    distinct() %>%
    dplyr::select(gene.name, dist, bait.id, oe.id, median.pi.diff, k.dmps)
}




# Promoters ---------------------------------------------------------------

extract_promoters_from_interactions <- function(interactions){

  # We'll keep the other end (enhancer) data sequestered
  oe.df <- data.frame(oeChr = seqnames(interactions),
                      oeStart = start(interactions),
                      oeEnd = end(interactions))

  tmp.df <- as.data.frame(interactions) %>%
    dplyr::select(-5:-1) %>%
    cbind(oe.df)

  makeGRangesFromDataFrame(tmp.df,
                           keep.extra.columns = T,
                           seqnames.field = "baitChr",
                           start.field = "baitStart",
                           end.field = "baitEnd")
}


# Big Lolly ---------------------------------------------------------------
# format_for_lolly <- function(data){
#   out <- data %>%
#     dplyr::transmute(chrom,
#                      start,
#                      end,
#                      name = ".",
#                      score = round(-log10(lfdr.from.ss)),
#                      strand = ".",
#                      thickStart = start,
#                      thickEnd = end,
#                      color = ifelse(pi.diff < 0, "0,119,154", "239,192,0"),
#                      lollySize = ifelse(lfdr.from.ss < 0.05, 5, 1))
#
#   # Make non-sig points grey
#   out[out$pValueLog < -log10(0.05), "color"] <- "220,220,220"
#   out
# }


# Export interactions ---------------------------------------------------------

export_significant_interactions_to_UCSC <- function(interactions.with.dmp){

  interactions.with.dmp %>%
    separate(oe.id, into = c("oeChr", "oeStart", "oeEnd")) %>%
    dplyr::transmute(sourceChrom = oeChr,
                     sourceStart = as.numeric(oeStart),
                     sourceEnd = as.numeric(oeEnd),
                     targetChrom = paste0("chr", baitChr),
                     targetStart = as.numeric(baitStart),
                     targetEnd = as.numeric(baitEnd)) %>%
    dplyr::mutate(chrom = sourceChrom,
                  chromStart = floor((sourceStart + targetStart) / 2),
                  chromEnd = floor((sourceEnd + targetEnd) / 2)) %>%
    dplyr::mutate(value = 5, exp = ".", color = "255,0,0",
                  name = ".", score = 500,
                  sourceName = "OtherEnd:enhancer", sourceStrand = ".",
                  targetName = "Bait:promoter", targetStrand = ".") %>%
    dplyr::select(all_of(ucsc.order.cols))
}


format_and_write_ucsc_interactions <- function(interactions.for.ucsc, file){

  out <-
    interactions.for.ucsc %>%
    drop_na() %>% # need distance to be non-missing to visualize...
    makeGRangesFromDataFrame(seqnames.field = "chrom",
                             start.field = "chromStart",
                             end.field = "chromEnd",
                             keep.extra.columns = T,
                             ignore.strand = T)

  # Extra processing
  genome(out) <- "hg38"

  lengths <- get_ucsc_seqlengths()

  # Allows us to get the indices to "duplicate" the
  # lengths the appropriate number of times. E.g., chr1 is 1mb,
  # then if there are five itneractions, we have 1mb, 1mb, 1mb, 1mb, 1mb
  ix <- as.vector(match(seqnames(out), names(lengths)))

  # Which interactions are in correct range
  keep.ix <- (end(out) <= lengths[ix]) & (start(out) <= lengths[ix])

  # Subset
  out.trimmed <- out[keep.ix, ]
  seqlengths(out.trimmed) <- lengths

  out.trimmed <- out.trimmed[,-1]
  
  # Write write write
  rtracklayer::export.bb(out.trimmed, file)
  file
}






get_summary_stats_pchic <- function(interactions.to.test, interactions.with.dmp){

  # N DMPs
  N.dmps.in.enhancer <- interactions.with.dmp[, 1:3] %>% distinct() %>% nrow()

  # Test set
  N.interactions.tested <- length(unique(interactions.to.test$interaction.id))
  N.promoter.baits.tested <- length(unique(interactions.to.test$bait.id))
  N.enhancer.oends.tested <- length(unique(interactions.to.test$oe.id))

  # Significant set
  N.interactions.sig <- length(unique(interactions.with.dmp$interaction.id))
  N.promoter.baits.sig <- length(unique(interactions.with.dmp$bait.id))
  N.enhancer.oends.sig <- length(unique(interactions.with.dmp$oe.id))

  return(list("N DMPs in at least one enhancer" = N.dmps.in.enhancer,
              "N interactions tested" = N.interactions.tested,
              "N promoter (baits) tested" = N.promoter.baits.tested,
              "N enhancer (other ends) tested" = N.enhancer.oends.tested,
              "N interactions significant" = N.interactions.sig,
              "N promoter (baits) significant" = N.promoter.baits.sig,
              "N enhancer (other ends) significant" = N.promoter.baits.sig))

}



# Integrate with RNAseq ---------------------------------------------------

unnest_interactions_by_gene <- function(inter){
  inter %>%
    as.data.frame() %>%
    mutate(gene.name = str_split(baitName, ";")) %>%
    unnest(gene.name)
}


get_common_genes_from_DM_interactions <- function(enhancers.with.dmp, promoters.with.dmp){
  ids <- intersect(enhancers.with.dmp$interaction.id, promoters.with.dmp$interaction.id)

  unnest_interactions_by_gene(enhancers.with.dmp) %>%
    dplyr::filter(interaction.id %in% ids) %>%
    pull(gene.name) %>%
    unique() %>%
    sort()
}


test_ranks_of_pchic_rnaseq <- function(diff.exp.data, dmps.in.enhancer.df){

  # Need one gene per row
  data.by.enhancers <- dmps.in.enhancer.df %>%
    dplyr::select(oe.id, baitName, median.pi.diff, k.dmps) %>%
    distinct() %>%
    unnest_interactions_by_gene() %>%
    dplyr::group_by(gene.name) %>%
    summarize(mean.dmps.per.enhancer = mean(k.dmps),
              N.enhancers = n_distinct(oe.id),
              median.pi.diff = median(median.pi.diff)) %>%
    distinct()


  # Joined by gene.name
  tmp <- inner_join(diff.exp.data, data.by.enhancers, by = "gene.name")

  xx <- abs(tmp$logFC)
  yy <- abs(tmp$median.pi.diff)

  return(list(
    "N enhancer genes" = length(unique(data.by.enhancers$gene.name)),
    "N common (DE, DME) genes" = length(unique(tmp$gene.name)),
    "Kendall's rank test p-value" = cor.test(xx, yy, method = "kendall")$p.val
  ))
}


compute_lambda_gc <- function(pp){
  round(median(qchisq(1 - pp, 1)) / qchisq(0.5, 1), 1)
}

# compute_lambda_robust <- function(x_vec){
#   (median(x_vec) / 0.675) ^ 2
# }


make_title_for_pvals_plots <- function(desc, lambda.gc){
  bquote(.(desc)~~"("*lambda[GC]==.(lambda.gc)*")")
}

plot_hist_from_pvals <- function(pp){

  lambda.gc <- compute_lambda_gc(pp)
  ti <- make_title_for_pvals_plots("Modeled p-values", lambda.gc)

  hist(pp,
       probability = T,
       breaks = 50,
       main = as.expression(ti),
       xlab = expression(p))
}


plot_qq_from_pvals <- function(pp){

  # Genomic in/de-flation
  lambda.gc <- compute_lambda_gc(pp)

  # Title
  ti <- make_title_for_pvals_plots("Raw p-values", lambda.gc)

  fastqq::qq(pp, main = as.expression(ti))
}


plot_pvals_routine <- function(pvals.data){
  pp.1 <- pvals.data$pval.Wald
  pp.2 <- pvals.data$pval

  # Run plotting (both coalculate lambda internally)

  # par(mfrow = c(1,2))
  plot_qq_from_pvals(pp.1)
  # plot_hist_from_pvals(pp.2)

}

run_and_save_pvals_plot <- function(pvals.data, file){
  pdf(file, width = 6, height = 6)
  plot_pvals_routine(pvals.data )
  dev.off()

  # Need to return a charcter vector
  file
}



# Coverage histograms -----------------------------------------------------

munge_sequencing_depth <- function(file, sample.ids){
  read_csv(file,
           show_col_types = F,
           col_names = c("sample_id", "Mean sequencing depth")) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    dplyr::filter(sample_id %in% sample.ids)
}


munge_effective_coverage <- function(file, sample.ids){

  df <- read_csv(file, show_col_types = F)

  # We'll need weighted averages based on number of positions per subject
  df$weight <- df$N.positions / sum(df$N.positions)

  # True sequencing coverage
  cov <- df %>%
    dplyr::select(-N.positions) %>%
    pivot_longer(-weight) %>%
    group_by(name) %>%
    summarize(Cov.bar = as.numeric(weight %*% value))

  names(cov) <-c("sample_id", "Mean effective coverage")
  cov %>%
    dplyr::filter(sample_id %in% sample.ids)
}



plot_coverage_hist_routine <- function(depth.df, cov.df, ofile){
  df <- full_join(cov.df, depth.df, by = "sample_id") %>%
    pivot_longer(-sample_id) %>%
    mutate(name = factor(name, levels = c("Mean sequencing depth", "Mean effective coverage")))

  # Plot this
  p <- df %>%
    ggplot(aes(x = value, y = after_stat(density))) +
    geom_histogram(bins = 20, color = "black", fill = "grey90") +
    theme_minimal() +
    facet_wrap(.~name, scales = "free") +
    ylab("Density") +
    xlab("") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "white", color = "white"))


 # cowplot::save_plot(ofile, p)
  pdf(ofile)
  print(p)
  dev.off()

}




# QC Stuff from JSONS -----------------------------------------------------

get_stats_from_json <- function(file){
  z <- fromJSON(file = file)

  sample_id = tools::file_path_sans_ext(basename(file))

  # Get the four values (total (L, R), unmapped (L, R))
  vv <- as.numeric(unlist(z$Reads))

  if(length(z$Reads) > 0){
    data.frame(sample_id = sample_id,
               TotalForward = vv[1],
               TotalReverse = vv[2],
               UnmappedForward = vv[3],
               UnmappedReverse = vv[4])
  } else {
    NA
  }
}

get_all_stats_from_dir <- function(dir){
  all.files <- list.files(dir, full.names = T, pattern = "*json")

  result <- lapply(X=all.files, FUN=get_stats_from_json)
  result <- result[!is.na(result)]

  do.call(rbind, result) %>%
    mutate(TotalReads = TotalForward + TotalReverse,
           TotalMappedReads = TotalReads - UnmappedForward - UnmappedReverse,
           MappingPercentage = TotalMappedReads / TotalReads)
}

# Write supplementary tables ----------------------------------------------


# Array comparisons -------------------------------------------------------

read_one_raw_M_Cov <- function(file){
  out <- fread(file) %>%
    dplyr::filter(coverage > 5) %>%
    dplyr::transmute(chrom,
                     start = chromStart,
                     P = methylated / coverage,
                     sample)
  print(file)
  return(out)
}

join_raw_M_Cov_files <- function(dir, common.ids){
  all.files <- file.path(dir, paste0("chr1.", common.ids, ".bed"))
  data.table::rbindlist(lapply(X = all.files, FUN = read_one_raw_M_Cov)) %>%
    pivot_wider(id_cols = start, names_from = sample, values_from = P) %>%
    dplyr::mutate(chr = "chr1", end = start + 2)
}

chunker_wrapper <- function(file){
  chunker(file,
          sep = "\t",
          chunksize = 500000,
          has_colnames = T,
          has_rownames = F)
}


make_P_from_M_Cov <- function(M, Cov, chr){
  P <- M[, -1] / Cov[ ,-1]

  # Make this explicit
  P[Cov[ ,-1] == 0] <- NA

  P$chr <- chr
  P$start <- M$chromStart
  P$end <- M$chromStart + 2

  makeGRangesFromDataFrame(P,
                           keep.extra.columns = T,
                           starts.in.df.are.0based = T)
}

read_and_merge_M_Cov_files <- function(file.prefix){
  chr <- basename(file.prefix)

  M.file <- paste0(file.prefix, ".M.bed")
  Cov.file <- paste0(file.prefix, ".Cov.bed")

  M_chunker <- chunker_wrapper(M.file)
  Cov_chunker <- chunker_wrapper(Cov.file)


  P.all.grs <- list()
  ix <- 1

  while(next_chunk(M_chunker) & next_chunk(Cov_chunker)){

    M <- get_table(M_chunker)
    Cov <- get_table(Cov_chunker)

    P.all.grs[[ix]] <- make_P_from_M_Cov(M, Cov, chr)
    ix <- ix + 1
  }

  return(do.call(c, P.all.grs))

}


read_850k_array <- function(file){
  array <- fread(file)

  array.gr <- makeGRangesFromDataFrame(array,
                                       keep.extra.columns = T,
                                       starts.in.df.are.0based = F)

  array.gr
}


combine_850k_array_with_wgms <- function(array.gr, P.gr){

  # Get common columns
  cc <- intersect(names(mcols(P.gr)), names(mcols(array.gr)))

  df.1 <- subsetByOverlaps(P.gr, array.gr)[ ,cc] %>%
    as.data.frame() %>%
    dplyr::mutate(tech = "WGMS")
  df.2 <- subsetByOverlaps(array.gr, P.gr)[ ,cc] %>%
    as.data.frame() %>%
    dplyr::mutate(tech = "EPIC")

  out.df <- rbind(df.1, df.2) %>%
    dplyr::select(-c(end, width, strand)) %>%
    pivot_longer(starts_with("X"), names_to = "sample", values_to = "methylation", ) %>%
    pivot_wider(values_from = "methylation", names_from = "tech", values_fn = list) %>%
    unnest(cols = c(WGMS, EPIC)) %>%
    distinct() %>%
    dplyr::mutate(sample = as.numeric(str_remove(sample, "X")))

  return(out.df)
}


merge_array_and_wgms_subroutine <- function(array.gr, file.prefix){
  P.gr <- read_and_merge_M_Cov_files(file.prefix)

  # We sample because otherwise it's impossible to
  # compute a correlation on 850k * 84 entries
  combine_850k_array_with_wgms(array.gr, P.gr)
}


merge_array_and_wgms_routine <- function(array.gr, dir){
  all.file.prefixes <- file.path(dir, paste0("chr", 1:22))

  set.seed(919)

  out <- lapply(X = all.file.prefixes,
                FUN = merge_array_and_wgms_subroutine,
                array.gr = array.gr) %>% do.call(rbind)

}


remove_dmps_from_merged <- function(merged, dmps.gr){

  # merged is a data frame with columns seqnames, start, WGMS, EPIC
  # dmps.gr is a genomic ranges


  zz <- makeGRangesFromDataFrame(
    dplyr::mutate(merged, end = start + 2),
    keep.extra.columns = TRUE
    )

  subsetByOverlaps(zz, dmps.gr, invert=TRUE, minoverlap = 2) %>%
    as.data.frame()

}


count_number_of_loci_dropped <- function(with_dmps, wo_dmps){
  tmp <- with_dmps %>%
    group_by(seqnames) %>%
    summarize(N = mean(parameter))
  N0 <- sum(tmp$N)

  tmp <- wo_dmps %>%
    group_by(seqnames) %>%
    summarize(N = mean(parameter))
  N1 <- sum(tmp$N)

  return(N1 - N0)

}

compute_wgms_vs_array_summary_stats <- function(chr_merged){

  res1 <- chr_merged %>%
    group_by(seqnames, sample) %>%
    do(broom::tidy(cor.test(.$WGMS, .$EPIC, use = "pairwise")))

  # res2 <- chr_merged %>%
  #   group_by(seqnames, sample) %>%
  #   do(broom::tidy(cor.test(z$WGMS, z$EPIC, use = "pairwise")))

  return(res1)
  # return(rbind(res1, res2))
}


plot_wgms_vs_array_hexbin <- function(wgms_vs_array.df){
  wgms_vs_array.df %>%
    drop_na() %>%
    sample_n(100000) %>%
    ggplot(aes(x = EPIC, y = WGMS)) +
    geom_hex(binwidth = c(0.01, 0.01), aes(fill = log(..density..))) +
    scale_fill_viridis(option = "plasma") +
    theme_classic()
}




# Combine combine combine -------------------------------------------------

curate_genes_by_dm_status <- function(dm.genes,
                                      dm.promoters.genes,
                                      dm.enhancers.genes,
                                      dex.genes,
                                      gwas.genes,
                                      array.genes){

  uniq.genes <- union(
      union(union(dm.genes, dex.genes), array.genes),
      union(union(dm.promoters.genes, dm.enhancers.genes), gwas.genes)
    )

  out <- data.frame(Gene_Symbol = uniq.genes) %>%
    dplyr::mutate(DiffMeth_Gene_Body = ifelse(Gene_Symbol %in% dm.genes, 1, 0),
                  DiffMeth_Promoter = ifelse(Gene_Symbol %in% dm.promoters.genes, 1, 0),
                  DiffMeth_Enhancer = ifelse(Gene_Symbol %in% dm.enhancers.genes, 1, 0),
                  DiffExp_Gene = ifelse(Gene_Symbol %in% dex.genes, 1, 0),
                  NatGen_GWAS_Gene = ifelse(Gene_Symbol %in% gwas.genes, 1, 0),
                  DiffMeth_EPICArray = ifelse(Gene_Symbol %in% array.genes, 1, 0)) %>%
    dplyr::mutate(TotalInclusions = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>%
    dplyr::arrange(-TotalInclusions)

  return(out)

}



get_paper_stats <- function(gene.body.enrichment, gene.lfdr){
  N <- sum(gene.body.enrichment$N.DMPs > 0 & gene.body.enrichment$lfdr < gene.lfdr)

  list("N genes with >0 DMPs + gene-wide lFDR < cutoff: " = N)
}



# General reading/casting functions -----------------------------------------------

get_natgen_genes <- function(ifile, exclude_APP = T, exclude_IGH = T){
  zz <- read_table(ifile, col_names = FALSE)$X1

  if (exclude_APP){
    zz <- zz[!(zz == "APP")]
  }

  if (exclude_IGH){
    zz <- zz[!stringr::str_detect(zz, "IGH")]
  }
  zz

}


my_write_csv <- function(data, file){
  write_csv(data, file = file)
  return(file)
}


get_pvals_data <- function(ifile){
  # READ pvals from
  data.table::fread(ifile, verbose = F) %>%
    dplyr::mutate(lfdr = lfdr.from.ss,
                  pval = p.from.ss,
                  pval.Wald = p.from.DSS,
		              pi.diff = pi.diff.mci.ctrl, ##### this is what you need to change for different datasets
                  y = -log10(lfdr)) %>%
    dplyr::select(-c(p.from.ss, p.from.zz, p.from.DSS,
                     lfdr.from.ss, lfdr.from.zz))
}

read_and_cast_madrid_data <- function(file){
  read_csv(file, show_col_types = F, skip = 1) %>%
    separate(hg38_coordinates, into = c("chr", "start")) %>%
    dplyr::mutate(start = as.numeric(start)) %>%
    dplyr::mutate(start = start - 1,
                  end = start + 2,
                  strand = Strand) %>%
    drop_na() %>%
    dplyr::transmute(chr, start, end,
                     gene.symbols = `Gene Symbol(s)`,
                     effect.size = -1 * `Mean Beta Difference (Control-AD)`) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T,
                             ignore.strand = T,
                             seqnames.field = "chr")
}


read_roubroeks_data <- function(ifile){
  readxl::read_xlsx(ifile, sheet = 1, skip = 5, col_names = F) %>%
    dplyr::select(c(1, 2, 7, 8)) %>%
    dplyr::rename("cpg" = 1, "locus" = 2, "difference" = 3, "tukey.p" = 4) %>%
    tidyr::separate(locus, into = c("chrom", "start"), sep = ": ") %>%
    dplyr::mutate(chrom = str_remove(chrom, " "),
                  start = as.numeric(start),
                  end = start + 2) %>%
    dplyr::filter(tukey.p < 0.01) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
}


lift_manifest <- function(manifest, chain){
  # Load the manifest file and get it in the right format
  anno <- as.data.frame(getAnnotation(manifest))

  # Filter and add start/end
  zz <- anno %>%
    dplyr::mutate(start = pos-1, end = pos+1) %>%
    makeGRangesFromDataFrame(starts.in.df.are.0based = T)

  # Finally, get to the right coordinate system
  unlist(rtracklayer::liftOver(zz, chain))
}


load_and_lift_450k <- function(chain){
  manifest <- IlluminaHumanMethylation450kanno.ilmn12.hg19
  lift_manifest(manifest, chain)
}


load_and_lift_EPIC <- function(chain){
  manifest <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19
  lift_manifest(manifest,chain)
}


filter_450k_annotation <- function(valid.cpgs){


}





# Constants ---------------------------------------------------------------


GENIC.MAPPING <-
  data.frame(
    Annotation = c("1to5kb", "3UTRs", "5UTRs", "exons", "introns", "promoters", "Other"),
    Annotation.clean = c("1-5 kb", "3' UTR", "5' UTR", "Exon", "Intron", "Promoter", "Other")
  )


CPG.MAPPING <-
  data.frame(
    Annotation = c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores", "Other"),
    Annotation.clean = c("Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")
  )

# Put factors in corret order, refelcting biology
GENIC.LEVELS <- c("1-5 kb", "Promoter", "5' UTR", "Exon", "Intron",  "3' UTR", "Other")
CPG.LEVELS <- c("Open sea", "CpG island", "CpG shelf", "CpG shore", "Other")



# Helper functions --------------------------------------------------------


annotate_loci_to_genic_parts <- function(data.gr, annotate.to = 'hg38_basicgenes'){

  if ("hg38_basicgenes" %in% annotate.to){
    MAPPING <- GENIC.MAPPING
    LEVELS <- GENIC.LEVELS
  } else {
    MAPPING <- CPG.MAPPING
    LEVELS <- CPG.LEVELS
  }

  print(MAPPING)

  # hg38_basicgenes or hg38_cpg
  my.anno <- annotatr::build_annotations(genome = 'hg38', annotations = annotate.to)

  # Run the annotation function
  data.annotated <- annotatr::annotate_regions(data.gr, my.anno, quiet = T)

  # Post-process after linking a CpG to a genomic feature
  anno.df <- data.frame(data.annotated) %>%
    dplyr::transmute(cpg.locus = paste0(seqnames, "-", start), annot.type) %>% # unique id for each CpG
    distinct() %>%  # get rid of instances where one CpG resides in multiple introns
    dplyr::mutate(Annotation = str_remove(annot.type, "hg38_genes_"))

  # Count how many annotate to none of these
  none.of.the.above <- subsetByOverlaps(data.gr,
                                        data.annotated,
                                        invert = T)

  # How many went to none?
  N.other <- length(none.of.the.above)

  # Add an extra row to anno.df; annotatr does not (as best as I know)
  # natively support counting "nones"
  other.df <- data.frame(cpg.locus = rep("none", N.other),
                         annot.type = rep("none", N.other),
                         Annotation = rep("Other", N.other))


  # Get the counts
  anno.counts.df <- rbind(anno.df, other.df) %>%
    group_by(Annotation) %>%
    dplyr::summarize(value = round(100 * n() / length(data.gr), 1))

  # Clean some strings and put in correct order
  anno.counts.cleaned.df <- anno.counts.df %>%
    left_join(MAPPING, by = "Annotation") %>%
    dplyr::mutate(Annotation.clean = factor(Annotation.clean, levels = LEVELS)) %>%
    arrange(Annotation.clean)

  # Add a string that looks like 'Promoter (10 %)'
  anno.counts.cleaned.df %>%
    dplyr::mutate(label = paste0(Annotation.clean,  " (", value, "%)")) %>%
    column_to_rownames("label")
}



cpg_annotation_routine <- function(data.gr){
  annotate_loci_to_genic_parts(
    data.gr,
    annotate.to = c("hg38_cpg_inter", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_shores")
  )
}




# Plotting ----------------------------------------------------------------


my_color <- 'd3.scaleOrdinal() .domain(["base", "not.genic", "genic", "other"])
.range(["#444444", "#800000", "#F2BA49", "#FFFF9F"])'


plot_sankey <- function(data, N.dmps){
  # Need rownames and column called "value"

  # Add a name for "Starting DMPs" (sink)
  ti <- htmltools::HTML(paste0("<strong><center>N = ", format(N.dmps, big.mark=","), " DMPs</center></strong>"))
  sink.name <- ""
  nodes <- data.frame(name = c(sink.name, rownames(data)))
  nodes$group <- as.factor(c("base", "not.genic", "not.genic", "not.genic",
                             "genic", "genic", "not.genic", "other"))


  # Munging for Sankey
  data$source <- 0
  data$target <- 1:nrow(data)


  # Return the object
  sankey.plot <- sankeyNetwork(Links = data,
                               Nodes = nodes,
                               Source = "source", Target = "target",
                               Value = "value", NodeID = "name",
                               fontSize = 12, nodeWidth = 20,
                               colourScale=my_color, NodeGroup="group")

  sankey.plot <- htmlwidgets::prependContent(sankey.plot, ti)

  sankey.plot
}

screenshot_sankey <- function(sankey.plot, filename){

  tmp.file <- ("tmp.sankey.html")
  saveNetwork(sankey.plot, file = tmp.file, selfcontained = F)

  # zoom increases the resolution it seems
  webshot::webshot(tmp.file, filename, vwidth = 500, vheight = 300, zoom=5)
  file.remove(tmp.file)

  return(filename)
}




# Volcano plotting --------------------------------------------------------


thin_large_pvals <- function(data){
  # Drop most of the non-significant points
  # For thinning non-sig points
  set.seed(919)

  ix <- sample(which(data$y < -log10(0.2)))
  ix.to.ignore <- ix[1:floor(length(ix) * 0.98)]

  # Subset
  data[-ix.to.ignore, ]
}

add_hyper_hypo_designation <- function(data, lfdr.cut=0.05){
  # Handle colors
  data$color <- "nonsig"
  data$color[data$lfdr < lfdr.cut & data$pi.diff > 0] <- "hyper"
  data$color[data$lfdr < lfdr.cut & data$pi.diff < 0] <- "hypo"

  # Order of color and order of pallete
  data$color <- factor(data$color, levels = c("nonsig", "hyper", "hypo"))
  data
}

plot_volcano <- function(data){
  my.colors <- get_hyper_hypo_colors()
  my.pal <- c("grey", my.colors$hyper, my.colors$hypo)

  data %>%
    drop_na() %>%
    ggplot(aes(x = pi.diff, y = y, color = color)) +
    geom_point(size = 0.3) +
    xlab("AD methylation % minus no-AD methylation %") +
    ylab(expression(-log[10](lFDR))) +
    scale_x_continuous(breaks = round(seq(-0.15, 0.15, 0.05), 2), lim = c(-0.15, 0.15)) +
    scale_color_manual(values = my.pal) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), alpha = 0.8) +
    theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.1),
          panel.grid.major.x = element_line(color = "grey", linewidth = 0.1),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 0),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"))
}


volcano_routine <- function(data){
  data %>%
    thin_large_pvals() %>%
    add_hyper_hypo_designation() %>%
    plot_volcano()
}



# Pi chart ----------------------------------------------------------------

dmp_pi_chart_routine <- function(data, file){
  my.colors <- get_hyper_hypo_colors()

  # How many DMPs
  N <- nrow(data)

  # Munge
  tally.df <- data %>%
#    dplyr::mutate(direction = ifelse(pi.diff.mci.ctrl < 0, "Hypomethylation", "Hypermethylation")) %>%
     dplyr::mutate(direction = ifelse(pi.diff < 0, "Hypomethylation", "Hypermethylation")) %>%
   group_by(direction) %>%
    summarize(prop = round(100 * n() / N)) %>%
    arrange(-prop) %>%
    dplyr::mutate(ypos = cumsum(prop) - 0.5*prop) %>%
    dplyr::mutate(label = paste0(direction, " (", prop, "%)"))

  p.pi <- ggplot(tally.df, aes(x="", y=prop, fill=label)) +
    geom_bar(stat="identity", width=1.5, color="white") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = "white"),
          plot.title = element_text(hjust = 0.5, size=20),
          legend.text=element_text(size=16),
          legend.position = "bottom",
          plot.margin = margin(t=5, r=5,b=5,l=5)) +
#    scale_fill_manual(name=NULL, values = c(my.colors$hyper, my.colors$hypo)) +
    scale_fill_manual(name=NULL, values = c("firebrick3", "mediumblue")) +
    guides(fill=guide_legend(ncol=1))
  pdf(file)
#  p.pi
  print(p.pi)
  dev.off()
}


# Plot CpG island, shore, shelf designation -------------------------------

plot_cpg_pi_chart <- function(data, file){
pdf(file)
  z <- data %>%
    dplyr::filter(value > 0.1) %>%
    rownames_to_column("Location") %>%
    ggplot(aes(fill = Location, x = "", y = value)) +
    geom_bar(width=1, stat="identity", color="white") +
    scale_fill_manual(name = NULL, values = pal_nejm()(4)) +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = "white"),
          plot.title = element_text(hjust = 0.5, size=20),
          legend.text=element_text(size=16),
          legend.position = "bottom",
          plot.margin = margin(t=5, r=5,b=5,l=5)) +
    guides(fill=guide_legend(ncol=2, reverse = F))
#    guides(fill=c("firebrick3","mediumblue"))
  print(z)
    dev.off()
}



# Stitch ------------------------------------------------------------------

stitch_birdseye_fig <- function(volcano.file, sankey.file, pi.file, cpg.file, out.file){

  volcano <- ggdraw() + draw_image(volcano.file)
  sankey <- ggdraw() + draw_image(sankey.file)
  pi <- ggdraw() + draw_image(pi.file)
  cpg <- ggdraw() + draw_image(cpg.file)

  topright <- plot_grid(pi, NULL, cpg, rel_widths = c(0.5, 0, 0.5),
                        labels = c("B","", "C"), nrow = 1)

  right <- plot_grid(topright, NULL, sankey,
                     ncol = 1,
                     rel_heights = c(0.5, 0, 0.5),
                     labels = c(NA, NA, "D"))

  z <- plot_grid(volcano, right,
            ncol = 2,
            rel_widths = c(0.5, 0.5),
            labels = c("A", NA)) +
    theme(plot.background = element_rect(fill = "white", color = "white")) +
    panel_border(remove = TRUE)

  cowplot::save_plot(z, filename = out.file)
  out.file

}



plot_tech_comp_scatter <-  function(merged){
  # Small helper to plot a scatter plot
  merged %>%
    sample_n(10000) %>%
    ggplot(aes(x = EPIC, y = WGMS)) +
    geom_point(alpha = 0.1) +
    theme_minimal() %>%
    return()
}

#

read_GODMC_meQTLs <- function(ifile){
  fread(ifile, select = c("cpg","snp","beta_a1","pval","cistrans","pval_are","pval_mre")) %>%
    dplyr::group_by(cpg) %>%
    dplyr::slice(1)
}





filter_450k_annotation <- function(valid.cpgs){

  # Load the manifest file and get it in the right format
  anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

  # Filter and add start/end
  anno.sig <- anno %>%
    dplyr::filter(Name %in% valid.cpgs) %>%
    dplyr::mutate(start = pos-1, end = pos+1)

  anno.sig.gr <- makeGRangesFromDataFrame(anno.sig, starts.in.df.are.0based = T)
  return(anno.sig.gr)
}


make_godmc_dmp_table <- function(
    dmp.and.mqtl,
    array.450.dmps,
    godmc.hg38.gr,
    array.450.hg38
    ){

  # Assume all 450k sites are included in universe
  N.universe <- length(array.450.hg38)

  x <- matrix(c(-1, -1, -1, -1), nrow = 2, ncol = 2)
  rownames(x) <- c("godmc_yes", "godmc_no")
  colnames(x) <- c("dmp_yes", "dmp_no")

  # Yes to both--DMP and meQTL
  x[1, 1] <- length(dmp.and.mqtl)

  # Yes to DMP, no to meQTL
  x[1, 2] <- length(array.450.dmps) - x[1, 1]

  # No to DMP, yes to meQTL
  x[2, 1] <- length(godmc.hg38.gr) - x[1, 1]

  # Inclusion-exclusion
  x[2, 2] <- length(array.450.hg38) - x[1, 2] - x[2, 1] + x[1, 1]

  tb <- as.table(x)
  tb
}

test_significance_of_meQTL_overlap <- function(N.meqtl.dmp, n.dmp, n.meqtl){
  # n.meqtl.dmp is the number of sites that are both DMPs and

  # See https://www.sciencedirect.com/science/article/pii/S0888754311001807?via%3Dihub
  # Bibikova 2014
  N.array <- 482421
  N.meqtl.dmp <- 10

}

compute_average_correlation <- function(df){
  (df$parameter %*% df$estimate) / sum(df$parameter)
}


count_n_dmps_multi_map_to_gene <- function(dmps, genes){
  seqlevelsStyle(dmps) <- "NCBI"
  seqlevelsStyle(genes) <- "NCBI"

  ov <- findOverlaps(dmps, genes)

  # Query is dmps in this case
  #check with max(queryHits(ov))
  data.frame(dmp.ix = queryHits(ov)) %>%
    dplyr::group_by(dmp.ix) %>%
    dplyr::summarize(N = n()) %>%
    return()
}

# tables_for_publication.R
# We include 5-6 tables in manuscript



process_and_write_dmps <- function(dmps.gr, desc, file){

  # Columns look like
  # Chrom, Start, End, AD - No AD, lFDR
  # chr1, 10, 12, 0.01, 0.02
  # ...

  data <- dmps.gr %>%
    as.data.frame() %>%
    dplyr::transmute(
      Chromosome = seqnames,
      "Start Coordinate (hg38)" = start,
      `End Coordinate (hg38)` = end,
#      `Estimated Methylation Difference (MCI minus CU)` = pi.diff.mci.ctrl,
       `Estimated Methylation Difference (MCI minus CU)` = pi.diff,
     `Local False-Discovery Rate (lFDR)` = lfdr)

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

#  return(file)
}


process_and_write_DM_genes <- function(gene.body.enrichment, desc, file){

  # Columns look like
  # Gene symbol, # CpGs, # DMPs, Harmonic Mean P-value, lFDR
  # BRCA1, 20, 1, 0.7, 0.99

  data <- gene.body.enrichment %>%
    transmute(
      `Gene Symbol` = gene_name,
      `Total Number of CpGs` = N.CpGs,
      `Number of Differentially Methylated Positions (DMPs)` = N.DMPs,
      `Harmonic Mean p-value` = HarmonicMeanPval,
      `Local False-Discovery Rate (lFDR)` = lfdr,
    )

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}



process_and_write_gene_ontology_terms <- function(dm.genes.go.df, desc, file){

  # First convert ENSEMBL IDs to genes...
  data <- dm.genes.go.df %>%
    transmute(`Gene Ontology (GO) Domain` = Ontology,
              `GO Term ID` = ID,
              `GO Description` = Description,
              `Gene Symbols` = GeneSymbols,
              `Gene Ratio` = GeneRatio,
              `BG Ratio` = BgRatio,
              `Local False-Discovery Rate (lFDR)` = p.adjust)

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)
}


format_and_write_dm_enhancers <- function(interactions.summary, de.genes, desc, file){

  data <- unnest_interactions_by_gene(interactions.summary) %>%
    dplyr::arrange(desc(k.dmps), gene.name) %>%
    transmute(
      `Number of Differentially Methylated Positions (DMPs)` = k.dmps,
      `Gene Symbol` = gene.name,
      `Promoter Locus (hg38)` = paste0("chr",baitChr, ":", baitStart, "-", baitEnd),
      `Enhancer Locus (hg38)` = paste0(oe.id),
      `Differentially Expressed` = ifelse(`Gene Symbol` %in% de.genes, "Differentially expressed", "Not differentially expressed")) %>%
    distinct() %>%
    arrange()


  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)
  return(file)
}


format_and_write_dm_promoters <- function(promoters.summary, de.genes, desc, file){

  data <- unnest_interactions_by_gene(promoters.summary) %>%
    dplyr::arrange(desc(k.dmps), gene.name) %>%
    transmute(
      `Number of Differentially Methylated Positions (DMPs)` = k.dmps,
      `Gene Symbol` = gene.name,
      `Promoter Locus (hg38)` = bait.id,
      `Enhancer Locus (hg38)` = paste0("chr", oeChr, ":", oeStart, "-", oeEnd),
      `Differentially Expressed` = ifelse(`Gene Symbol` %in% de.genes, "Differentially expressed", "Not differentially expressed")) %>%
    distinct()


  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)
  return(file)
}



process_and_write_nature_genetics_list <- function(NG.tally.df,
                                                   natgen.genes,
                                                   genes.w.dm.enhancers,
                                                   genes.w.dm.promoters, desc, file){


  # Need to expand current dataset
  genes.zero.dmps <- setdiff(natgen.genes, NG.tally.df$symbol)
  tmp <- data.frame(genes.zero.dmps, 0)
  colnames(tmp) <- colnames(NG.tally.df)

  # COmbine non-zero with zero counts
  data <- rbind(NG.tally.df, tmp) %>%
    dplyr::filter(!str_detect(symbol, "IGH")) %>%
    arrange(-N.DMPs) %>%
    dplyr::mutate(
      dm.enhancer = ifelse(symbol %in% genes.w.dm.enhancers, "1+ DMPs in 1+ enhancers", "None"),
      dm.promoter = ifelse(symbol %in% genes.w.dm.promoters, "1+ DMPs in 1+ promoters", "None"))
  head(data)

  colnames(data) <- c("AD Genomic Risk Locus Symbol",
                      "Number Of Differentially Methylated Positions (DMPs) Within 25kb",
                      "Differential Methylation in Enhancer",
                      "Differential Methylation in Promoter")

  wb <- create_wb_with_description(desc)
  wb <- append_data_to_workbook(wb, data)

  openxlsx::saveWorkbook(wb, file, overwrite = T)

  return(file)

}


theme_meAD <- function(){
  font <- "Arial"   #assign font family up front

  theme_minimal(base_size = 18) %+replace%    #replace elements we want to change

    theme(

      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks

      #since theme_minimal() already strips axis lines,
      #we don't need to do that again

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      #
      # plot.caption = element_text(           #caption
      #   family = font,            #font family
      #   size = 9,                 #font size
      #   hjust = 1),               #right align
      #
      # axis.title = element_text(             #axis titles
      #   family = font,            #font family
      #   size = 10),               #font size
      #
      # axis.text = element_text(              #axis text
      #   family = font,            #axis famuly
      #   size = 9),                #font size
      #
      # axis.text.x = element_text(            #margin for axis text
      #   margin=margin(5, b = 10))

    )
}

### End
