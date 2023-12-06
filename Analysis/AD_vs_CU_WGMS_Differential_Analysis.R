# REDO LOAD v CONTROL (using all samples for fitting)

### C.E.B. wrote the vast majority of these scripts
### Simply modified them to loop through each chromosome

# 0_prepare_DSS_chrom.R
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
for (i in 1:length(chr)) {


parser <- ArgumentParser()
parser$add_argument("--idir", default= "/media/Data/WGBS/LOAD_MCI/MCovRaw/", help='Directory to run DSS on')
parser$add_argument("--odir", default= "/media/Data/WGBS/LOAD_MCI/Inputs/CONTROL_MCI_LOAD_New/", help='Directory to run DSS on')
parser$add_argument("--blood_file", default="/media/Data/WGBS/LOAD_MCI/bloodCellPropsDXM.csv")
parser$add_argument("--master_file", default="/media/Data/WGBS/LOAD_MCI/masterSamplesheetFilteredSamples.csv")
#parser$add_argument("--chr", default= "chr22", help='Chromosome to run on')
parser$add_argument("--chr", default= chr[[i]], help='Chromosome to run on')
parser$add_argument("--median_threshold", default=5, help="Keep sites with at least () median coverage")
parser$add_argument("--percent_nonzero_threshold", default=0.5, help="Keep sites with at least () minimum coverage")
parser$add_argument("--top_percent_most_variable", default=0.05, help="Keep top (5%) most variable in PC computation")
parser$add_argument("--save_pivots", action="store_true", help= "Save pivoted M and Cov matrices? Usually used when you want to check.")
args <- parser$parse_args()

############ FILES ##############
odir <- args$odir
pc.odir <- file.path(odir, "PrinComps/")

# Two birds one stone
dir.create(pc.odir, showWarn = F, recursive = T)

ofile <- file.path(odir, paste0("filtered-DSS-inputs-", args$chr, ".RData"))
if (file.exists(ofile)) {
    stop(paste0(ofile, " already exists; delete if you want to run 0-prepare-DSS-inputs.R"))
}


#TODO--send this code to generate master... don't want to have to repeat this
master.df <- read.csv(args$master_file) %>%
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::select(-starts_with("PC")) %>%
    dplyr::mutate(
        diagnostic_group_coded = 
	case_when(diagnostic_group == "CONTROL" ~ 0,
	diagnostic_group == "MCI" ~ 1,
	diagnostic_group == "LOAD" ~ 2))

# Get and remove control samples that aren’t CU anymore
#removeIDs <- read.table("/media/Data/WGBS/LOAD_MCI/removeSamples.txt",header=T)
#master.df <- subset(master.df,!(study_id %in% removeIDs$sampleID))

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
Cov <- pivot_me(data, "coverage")

# Set to zero so we can compute summary stats
Cov.zeroed <- data.table::copy(Cov) 
Cov.zeroed[is.na(Cov.zeroed)] <- 0

# min.coverage <- apply(Cov.zeroed, FUN=min, MARGIN=1)
percent.nonzero <- 1 - (rowSums(Cov.zeroed == 0) / ncol(Cov.zeroed))
med.coverage <- apply(Cov.zeroed, FUN=median, MARGIN=1)

keepix <- ((percent.nonzero >= args$percent_nonzero_threshold) & 
           (med.coverage >= args$median_threshold))

percent.keep <- round(100 * sum(keepix) / length(keepix), 2)
paste0("Keeping ", percent.keep, "%")

# Save the M and Cov matrices if requested
if (args$save_pivots){
    ofile <- file.path(odir, paste0("pivot-", args$chr, ".RData"))
    print(paste0("Saving ", ofile))

    save(M, Cov, percent.keep, file = ofile)
}

# No need to impute if there's a minimum filter
M.filt <- data.table::copy(M)[keepix, ]
Cov.filt <- data.table::copy(Cov)[keepix, ]

# Number of samples
N <- ncol(M.filt)

create_filler_mask <- function(DT, ctrl.samples, mci.samples, load.samples, N){

    # N is number of columns (samples)
    row.means.load <- round(rowMeans(DT[ , load.samples], na.rm = T))
    row.means.mci <- round(rowMeans(DT[ , mci.samples], na.rm = T))
    row.means.ctrl <- round(rowMeans(DT[ , ctrl.samples], na.rm = T))

    # Expand to make it a matrix we can mask
    filler <- DT
    filler[ , mci.samples] <- row.means.mci    
    filler[ , load.samples] <- row.means.load
    filler[ , ctrl.samples] <- row.means.ctrl

    filler
}

# LOAD/Control row means expanded to be the same shape as M 
M.filler <- create_filler_mask(M.filt, ctrl.samples, mci.samples, load.samples, N)
Cov.filler <- create_filler_mask(Cov.filt, ctrl.samples, mci.samples, load.samples, N)

# Masks for indexing
mask <- (Cov.filt == 0 | is.na(Cov.filt))
percent.of.kept.imputed <- round(100 * sum(mask) / (nrow(M.filt) * ncol(M.filt)), 2)

# Fill withe means
M.filt[mask] <- M.filler[mask]
Cov.filt[mask] <- Cov.filler[mask]
invisible(gc())

#########################################
###### Compute Principal Components #####
#########################################

P <- t(M.filt / Cov.filt)
dim(P)

P.vars <- apply(P, 2, stats::var)

# Get a numeric cut from the percentile
cut <- quantile(P.vars, 1 - args$top_percent_most_variable)

# Subet methylation matrix  
P.top <- P[ ,(P.vars >= cut)]
dim(P.top)

# Compute PCs
# P.top should be (and is) wide
pca.out <- gmodels::fast.prcomp(P.top, scale. = T, center = T)
PC.df <- as.data.frame(pca.out$x) 
PC.df$sample_id <- rownames(PC.df)


pc.ofile <- file.path(pc.odir, paste0("pcs-", args$chr, ".pdf"))
pdf(pc.ofile)
plot(PC.df$PC1, PC.df$PC2,
    xlab = "PC1", ylab = "PC2")
dev.off()


#########################################
########## Samplesheet munging ##########
#########################################

#Load data with a bit of munging. Column names are better standardized on the full dataset
blood.df <- read.csv(args$blood_file) %>% 
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = as.character(sample_id))

#Pull phenotypes, PCs, and blood into one design matrix
#Subset to only have sample_ids also contained in M/Cov
df <- dplyr::inner_join(blood.df, master.df, by = "sample_id") %>%
        dplyr::inner_join(PC.df, by = "sample_id")

#M(ethylated) reads and Cov(erage) from sequencing
# after processing
M <- M.filt
Cov <- Cov.filt

# Correct sample IDs in the correct order
ordered_sample_ids <- c(intersect(df$sample_id, names(M)))
# Don't get rid of the distinct
df <- df %>% 
    dplyr::filter(sample_id %in% ordered_sample_ids) %>% 
    dplyr::distinct()

# Rearrange columns to match data frame
M <- M %>% dplyr::select(all_of(ordered_sample_ids))
Cov <- Cov %>% dplyr::select(all_of(ordered_sample_ids))

#Check that ordering of samples in M/Cov is the same as in samplesheet
#should be in ascending order
if (all(df$sample_id == names(M))){
    print("Sample order is correct")
} else {
    warning("Sample order wrong!")
}

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M)), 
            pos = as.numeric(rownames(M)),
            M = as.matrix(M), 
            Cov = as.matrix(Cov), 
            sampleNames = names(M))

# Order positions out of abundance of caution
bs <- orderBSseq(bs)

# Not neccessary for running serial, but when parallel
# it helps to clean up
invisible(gc())

# Derived the output directory from the name of the input
save(bs, df, pca.out, percent.keep, percent.of.kept.imputed, med.coverage, file = ofile)
print(paste0("Wrote out ", ofile))
}}
loopByChromosome0(chrs)


######################################## 
######################################## 
######################################## 
######## PART 2: ELECTRIC BOOGALOO ########
######################################## 
######################################## 
######################################## 

# 01-find-DMRs.R
suppressPackageStartupMessages({
    library(argparse)
    library(DSS)
    library(magrittr)
    library(dplyr)
    library(stringr)
    library(parallel)
})

# loop by chromosome
# get chromosome list
chrs <- list.files("/media/Data/WGBS/LOAD_MCI/MCovRaw/BismarkFormatted/chromosomes")
chrs <- chrs[1:22]
loopByChromosome1 <- function(chr) {
for (i in 1:length(chr)) {

cat(chr[[i]])
cat("\n")
cat("Setting up arguments . . .\n")

parser <- ArgumentParser()
parser$add_argument("--idir", default= "/media/Data/WGBS/LOAD_MCI/Inputs/CONTROL_MCI_LOAD_New/")
parser$add_argument("--odir", default= "/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/LOAD_CONTROL/test-diagnostic-group-coded/")
#parser$add_argument("--chr", default= "chr22")
parser$add_argument("--chr", default= chr[[i]])
parser$add_argument("--baseline_covariates", default= "bmi,sex,age_at_visit,type1,type2,type3,type4,type5,PC1,PC2")
parser$add_argument("--test_covariate", default= "diagnostic_group")
parser$add_argument("--smoothing", default= 150, help= 'Width of smoothing window')
args <- parser$parse_args()

ifile <- file.path(args$idir, paste0("filtered-DSS-inputs-", args$chr, ".RData"))

###########################################################
##################### Directory Manipulation ##############
###########################################################

cat("Reading in data . . .\n")

load(ifile)

# Temp patch, fixed in 0-prepare
# Should not have been storing a `args`
args <- parser$parse_args()

# Put the two paths together and create the output
odir <- args$odir
dir.create(odir, showWarn=F, recursive=T)

file_name <- str_remove(str_replace(basename(ifile), "input", "output"), "filtered-")
outname <- file.path(odir, file_name)

###########################################################
##################### DATA PREP ###########################
###########################################################
cat("Setting up formula and matrix . . .\n")
covariates_split <- stringr::str_split(args$baseline_covariates, ",")[[1]]
dss.formula <- formula(paste0("~", args$test_covariate, "+", paste(covariates_split, collapse = "+")))

# Pull out just the covariates we'll use
design.df <- df %>% 
    tibble::column_to_rownames("sample_id") %>%
    dplyr::select(all_of(c(args$test_covariate, covariates_split))) %>%
    tidyr::drop_na()

# If there are NAs in any of the covariates specified, be sure to drop
#corresponding sample in bs object
ix <- colnames(bs) %in% rownames(design.df)
bs.sub <- bs[ , ix]

if (all(colnames(bs.sub) == rownames(design.df))){
    print("After filtering NAs, rownames in design match colnames in bs")
} else {
    warning("rownames in design != colnames in bs")
}

# Derive some parameters
smooth = TRUE
if (args$smoothing == 0){
    smooth = FALSE
}

# Fit models
cat("Fitting models . . .\n")

dml.fit <- DMLfit.multiFactor(bs.sub, design = design.df, smoothing = smooth, 
    smoothing.span = args$smoothing, formula = dss.formula)

cat("\n")

invisible(gc())

# Test covariate (to start this will be LOAD status--diagnostic_group_coded)
cat("Running differential analysis . . .\n")
test.var <- args$test_covariate
test.result <- DMLtest.multiFactor(dml.fit, coef = "diagnostic_groupLOAD")

#To measure effect size, other beta coefficients...
beta.df <- data.frame(dml.fit$fit$beta)
colnames(beta.df) <- colnames(dml.fit$X)

design.mat <- model.matrix(dss.formula, design.df)

cat("Saving data . . .\n")

save(list = intersect(ls(), c("beta.df", "test.result", "design.df", "dss.formula", "test.var", "design.df")), 
    file = outname)
}}
loopByChromosome1(chrs)



######################################## 
######################################## 
######################################## 
############### PART 3 ##################
######################################## 
######################################## 
######################################## 



# 2-extract-pi-and-pvals.R
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(argparse)
    library(DSS)
    library(fdrtool)
    library(GenomicRanges)
    library(Rcpp)
})

Rcpp::sourceCpp("/media/Data/WGBS/LOAD_MCI/matrix_multiply.cpp")

parser <- ArgumentParser()
parser$add_argument("--idir", default= "/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/LOAD_CONTROL/test-diagnostic-group-coded/", help='Where are the models stored')
parser$add_argument("--odir_pis", default= "/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/LOAD_CONTROL/Outputs/Pis/", help='Where are the models stored')
parser$add_argument("--odir_summaries", default= "/media/Data/WGBS/LOAD_MCI/Results/CONTROL_MCI_LOAD/LOAD_CONTROL/Outputs/Summaries/", help='Where are the models stored')
parser$add_argument("--genomic_control", action="store_true")
args <- parser$parse_args()

# Output directories
odir.pis <- args$odir_pis
dir.create(odir.pis, showWarn=F)
print(paste0("Will write pis ", odir.pis))

odir.summaries <- args$odir_summaries
dir.create(odir.summaries, showWarn=F)
print(paste0("Will write pvals.df ", odir.summaries))

# Files to eventually load
all.files <- list.files(args$idir, pattern = "*RData", full=T)

compute_pi_for_cpg_i <- function(i, X, beta.mat){
    # helper function that calls KMP's matrix multiply script
    #i is CpG index
    b <- beta.mat[i, ]

    # Inverse link function
    # y = arcsin(2*X\beta - 1)
    tmp <- (sin(multiply(as(X,'sparseMatrix'), as(b, 'sparseMatrix'))) + 1) / 2
    as.numeric(tmp)
}


compute_pis <- function(X, beta.mat, test.result){

    N <- nrow(beta.mat)

    out <- mclapply(FUN = function(i) {
        compute_pi_for_cpg_i(i, X, beta.mat)
        }, 
        1:N,mc.cores=12
    )

    pi.mat <- do.call(rbind, out)

    # Assign column names based on sample id
    colnames(pi.mat) <- rownames(X)

    # Add chr, start, end as leading columns
    class(test.result) <- "data.frame" # remove DSS class
    bed.cols <- dplyr::transmute(test.result, chr, start = pos, end = pos + 2)

    pi.df <- data.table(cbind(bed.cols, pi.mat), rownames = NULL)
    pi.df
}

compute_pi_diffs <- function(X, design.df, pi.df){
    
    load.samples <- rownames(X)[design.df$diagnostic_group == "LOAD"]
    mci.samples <- rownames(X)[design.df$diagnostic_group == "MCI"]
    control.samples <- rownames(X)[design.df$diagnostic_group == "CONTROL"]

    pi.diff.load.ctrl <- rowMeans(pi.df[ , ..load.samples]) - 
                rowMeans(pi.df[ , ..control.samples])
    pi.diff.load.mci <- rowMeans(pi.df[ , ..load.samples]) - 
                rowMeans(pi.df[ , ..mci.samples])
    pi.diff.mci.ctrl <- rowMeans(pi.df[ , ..mci.samples]) - 
                rowMeans(pi.df[ , ..control.samples])
     pi.diff <- cbind(pi.diff.load.mci,pi.diff.load.ctrl,pi.diff.mci.ctrl)

}
make_pi_oname <- function(ifile){
    z <- stringr::str_replace(
            stringr::str_replace(basename(ifile), "output", "pi"), 
        "RData", "bed")
    file.path(odir.pis, z)
}


run_pi_routine <- function(ifile){
    ofile <- make_pi_oname(ifile)

    if (!file.exists(ofile)){
        load(ifile)

        # While data is loaded, compute pis and 
        X <- model.matrix(dss.formula, design.df)
        beta.mat <- as(beta.df, "matrix")

        pi.df <- compute_pis(X, beta.mat, test.result)
        print(paste0("Writing out ", ofile))

        fwrite(pi.df, ofile, sep="\t", verbose =F)
        gc()

    } else {
        print(paste0(ofile, " already exists, skipping"))
    }
}

lapply(FUN=run_pi_routine, all.files)

# test.cohort has chr,pos,stat,pvals,fdrs as columns and is DSS-sepcific
load_data_with_pi_diff <- function(ifile){
    load(ifile)

    # Very inefficient to load back in but oh well
    pi.df <- fread(make_pi_oname(ifile))

    X <- model.matrix(dss.formula, design.df)
    pi.diffs <- compute_pi_diffs(X, design.df, pi.df)

    # covariates and p-values
    df <- cbind(test.result, beta.df,pi.diffs) 
    #df$pi.diff <- pi.diffs
   
    gc()

    return(df)
}
df <- do.call("rbind", lapply(all.files, load_data_with_pi_diff)) 

# Genomic inflation (test stats from DSS are N(0,1))
# https://github.com/haowulab/DSS/blob/55cb1a05738a2ed02717af50b7b52828bc6b508d/R/DML.multiFactor.R#L192
chisq <- df$stat^2
lambda <- median(chisq) / qchisq(0.5, 1)
stat.crct <- chisq / lambda

# Correct by adjusting with genomic inflation 
if (args$genomic_control){
    pp <- pchisq(stat.crct, df=1, lower.tail=FALSE)
} else {
    pp <- df$pvals
}

# Stats themselves
ss <- df$stat
# Z normalized stats
zz <- (ss - mean(ss)) / sd(ss)

##############################################################
################ FDR CACLULATION #############################
##############################################################
# False discovery rate, and 
pp.out <- fdrtool(pp, statistic = "pvalue", plot = FALSE)
ss.out <- fdrtool(ss, statistic = "normal", plot = FALSE)
zz.out <- fdrtool(zz, statistic = "normal", plot = FALSE)

# Rename columns, but don't sort yet
pvals.df <- df %>% 
    dplyr::mutate(chr = chr, start = pos , end = pos + 2) %>% 
    #dplyr::select("chr", "start", "end", 7, "pi.diff", "stat")
    dplyr::select("chr", "start", "end", 7, "pi.diff.load.mci", "pi.diff.load.ctrl", "pi.diff.mci.ctrl", "stat")

# Add in the pval-type stuff we want
pvals.df$p.from.ss <- ss.out$pval
pvals.df$p.from.zz <- zz.out$pval
pvals.df$p.from.DSS <- df$pvals

# lFDR
pvals.df$lfdr.from.ss <- ss.out$lfdr
pvals.df$lfdr.from.zz <- zz.out$lfdr
ofile <- file.path(odir.summaries, "pvals.bed")
fwrite(pvals.df, ofile, sep="\t")
