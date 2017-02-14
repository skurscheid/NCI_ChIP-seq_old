#!/usr/bin/env Rscript

# The R command line script simulates reads for a TF ChIP experiment using
# the R package ChIPsim by Peter Humburg
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 14/02/2017
# Last change: 14/02/2017


## ------------------------------------------------------------------------
## Start clock
t0 <- Sys.time();


## ------------------------------------------------------------------------
## Load libraries
suppressMessages(library(optparse));    # For Python-style command line args
suppressMessages(library(Biostrings));  # For DNAStringSet
suppressMessages(library(ChIPsim));     # Main package for ChIP-seq data simulation
suppressMessages(library(ggplot2));     # For plotting
suppressMessages(library(actuar));      # For pareto distribution
suppressMessages(library(reshape2));    # For melt (convert wide df to long df)


## ------------------------------------------------------------------------
## Parse command line arguments
option_list <- list(
    make_option(
        c("-r", "--ref"),
        type = "character",
        default = NULL,
        help = "Reference fasta file (fa. or .fa.gz)",
        metavar = "character"),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "ChIPsim",
        help = "Output name (without file extension)",
        metavar = "character"),
    make_option(
        c("--seed"),
        type = "integer",
        default = NULL,
        help = "Set fixed seed for random number generator [default %default]",
        metavar = "integer"),
    make_option(
        c("--bindLength"),
        type = "integer",
        default = 50,
        help = "Length of binding site in bp [default %default]",
        metavar = "integer"),
    make_option(
        c("--backLength"),
        type = "integer",
        default = 500,
        help = "Length of background site in bp [default %default]",
        metavar = "integer"),
    make_option(
        c("--meanFragLength"),
        type = "integer",
        default = 200,
        help = "Mean DNA fragment length in bp [default %default]",
        metavar = "integer"),
    make_option(
        c("--minFragLength"),
        type = "integer",
        default = 150,
        help = "Minimum DNA fragment length in bp [default %default]",
        metavar = "integer"),
    make_option(
        c("--maxFragLength"),
        type = "integer",
        default = 250,
        help = "Maximum DNA fragment length in bp [default %default]",
        metavar = "integer"),
    make_option(
        c("--readLength"),
        type = "integer",
        default = 75,
        help = "Read length in bp [default %default]",
        metavar = "integer"),
    make_option(
        c("-n", "--nReads"),
        type = "integer",
        default = 1e5,
        help = "Number of reads [default %default]",
        metavar = "integer")
    );
opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);
if (is.null(args$ref)) {
    print_help(opt_parser);
    stop("At least reference file must be supplied.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Timestamp function
ts <- function() {
    return(format(Sys.time(), "[%a %b %d %Y %H:%M:%S]"));
}


## ------------------------------------------------------------------------
## Global variables
genome <- args$ref;
seed <- args$seed;
backgroundLength <- args$backLength;
bindingLength <- args$bindLength;
outputName <- args$out;
meanFragmentLength <- args$meanFragLength;
minFragmentLength <- args$minFragLength;
maxFragmentLength <- args$maxFragLength;
readLength <- args$readLength;
nReads <- args$nReads;
cat(sprintf("%s Parameter summary\n", ts()));
cat(sprintf(" Reference          = %s\n", genome));
cat(sprintf(" outputName         = %s\n", outputName));
cat(sprintf(" seed               = %s\n", seed));
cat(sprintf(" bindingLength      = %i\n", bindingLength));
cat(sprintf(" backgroundLength   = %i\n", backgroundLength));
cat(sprintf(" meanFragmentLength = %i\n", meanFragmentLength));
cat(sprintf(" minFragmentLength  = %i\n", minFragmentLength));
cat(sprintf(" maxFragmentLength  = %i\n", maxFragmentLength));
cat(sprintf(" readLength         = %i\n", readLength));
cat(sprintf(" nReads             = %i\n", nReads));


## ------------------------------------------------------------------------
## Read in reference
genome <- readDNAStringSet(genome);
refLength <- width(genome);
cat(sprintf(
    "%s Read reference sequence \"%s\" of length %i\n",
    ts(),
    names(genome),
    refLength));
if (refLength < backgroundLength + bindingLength) {
    stop("Reference sequence must be longer than bindingLength + backgroundLength.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Define transition probabilities and initial state
transition <- list(
    Binding = c(Background = 1),
    Background = c(Binding = 0.05, Background = 0.95)
);
transition <- lapply(transition, "class<-", "StateDistribution");
init <- c(Binding = 0, Background = 1);
class(init) <- "StateDistribution";


## ------------------------------------------------------------------------
## Define background binding strength distribution
backgroundFeature <- function(
    start,
    length = backgroundLength,
    shape = 1,
    scale = 20) {
    weight <- rgamma(1, shape = 1, scale = 20);
    params <- list(start = start, length = length, weight = weight);
    class(params) <- c("Background", "SimulatedFeature");
    return(params);
}

## ------------------------------------------------------------------------
## Define TF binding strength distribution
bindingFeature <- function(
    start,
    length = bindingLength,
    shape = 1,
    scale = 20,
    enrichment = 5,
    r = 1.5) {
    stopifnot(r > 1);
    avgWeight <- shape * scale * enrichment;
    lowerBound <- (r - 1) * avgWeight;
    weight <- actuar::rpareto1(1, r, lowerBound);
    params <- list(start = start, length = length, weight = weight);
    class(params) <- c("Binding", "SimulatedFeature");
    return(params);
}


## ------------------------------------------------------------------------
## Generate feature sequences using Markov chain model
## Use fixed seed if seed from command line args is not null
if (!is.null(seed)) {
    cat(sprintf("%s Set fixed seed %i\n", ts(), seed));
    set.seed(seed);
}
cat(sprintf("%s Running MC model...\n", ts()));
generator <- list(
    Binding = bindingFeature,
    Background = backgroundFeature);
features <- ChIPsim::placeFeatures(
    generator,
    transition,
    init,
    start = 0,
    length = refLength,
    globals = list(shape = 1, scale = 20),
    experimentType = "TFExperiment",
    lastFeat = c(Binding = FALSE, Background = TRUE));


## ------------------------------------------------------------------------
## Stop if less <=2 binding sites
t <- table(sapply(features, class));
nBinding <- t[which(names(t) == "Binding")];
nBackground <- t[which(names(t) == "Background")];
cat(sprintf("%s Results of MC model\n", ts()));
cat(sprintf(" %10i binding sites\n", nBinding));
cat(sprintf(" %10i background sites\n", nBackground));
if (nBinding < 3) {
    stop("MC model gave <=2 binding sites. Decrease binding site lengths or change reference.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Plot the distribution of TF and background binding site strengths
df <- cbind.data.frame(
    weight = sapply(features, "[[", "weight"),
    type = sapply(features, function(x) class(x)[1]));
gg <- ggplot(subset(df, weight <= 500), aes(weight, fill = type));
gg <- gg + geom_density(alpha = 0.3);
gg <- gg + theme_bw();
gg <- gg + scale_fill_brewer(palette = "Dark2", name = "Binding site type");
gg <- gg + xlim(0, 500);
gg <- gg + labs(
    x = "Binding site strength",
    y = "Density",
    title = "Binding site strength distribution",
    subtitle = sprintf(
        "No of TF binding sites = %i, No. of background sites = %i",
        nBinding,
        nBackground));
ggsave(
    filename = sprintf("distr_binding_site_strength_%s.pdf", outputName),
    width = 8,
    height = 4,
    gg);


## ------------------------------------------------------------------------
## Define binding site densities
featureDensity.Binding <- function(feature, ...) {
  featDens <- numeric(feature$length);
  featDens[floor(feature$length / 2)] <- feature$weight;
#  featDens[floor(feature$length / 2)] <- feature$weight * feature$length;
  return(featDens);
}
featureDensity.Background <- function(feature, ...) {
  featDens <- numeric(feature$length);
#  featDens[] <- feature$weight;
  featDens[] <- feature$weight / bindingLength;
  return(featDens);
}


## ------------------------------------------------------------------------
## Convert features to densities
cat(sprintf("%s Converting features to site densities...\n", ts()));
dens <- ChIPsim::feat2dens(features, length = refLength);


## ------------------------------------------------------------------------
## Define DNA fragment size distribution
fragLength <- function(x, minLength, maxLength, meanLength, ...) {
  sd <- (maxLength - minLength)/4;
  prob <- dnorm(minLength:maxLength, mean = meanLength, sd = sd);
  prob <- prob/sum(prob);
  return(prob[x - minLength + 1]);
}


## ------------------------------------------------------------------------
## Convert binding site densities to DNA fragment densities
cat(sprintf("%s Converting site densities to DNA fragment densities...\n", ts()));
readDens <- ChIPsim::bindDens2readDens(
  dens,
  fragLength,
  bind = bindingLength,
  minLength = minFragmentLength,
  maxLength = maxFragmentLength,
  meanLength = meanFragmentLength);


## ------------------------------------------------------------------------
## Sample reads
cat(sprintf("%s Sampling reads...\n", ts()));
readLoc <- ChIPsim::sampleReads(readDens, nreads = 1e5);


## ------------------------------------------------------------------------
## Define uniform read quality function
randomQualityPhred33 <- function(read, ...) {
  # Character vector of symbols for the Phred+33 quality encoding scale
  rangePhred33 <- unlist(strsplit(rawToChar(as.raw(33:126)), ""));
  # Uniform-randomly sample qualities
  paste(sample(rangePhred33, length(read), replace = TRUE), collapse = "");
}

# We need to make sure that readLoc + readLen <= refLength for both strands
readLoc[[1]] <- readLoc[[1]][which(readLoc[[1]] + readLength <= refLength)];
readLoc[[2]] <- readLoc[[1]][which(readLoc[[1]] - readLength > 0)];

# Create names
nreads <- sapply(readLoc, length);
names <- list(fwd = sprintf("read_%s_fwd_%s", outputName, seq(nreads[1])),
              rev = sprintf("read_%s_rev_%s", outputName, seq(nreads[2])));

# Write to FASTQ
# Uncomment for output
cat(sprintf("%s Writing reads to fastq file...\n", ts()));
pos2fastq(readLoc,
          names = names,
          sequence = genome[[1]],
          qualityFun = randomQualityPhred33,
          errorFun = readError,
          readLen = readLength,
          file = sprintf("simul_%s.fastq", outputName));


## ------------------------------------------------------------------------
## Done
cat(sprintf("%s Done (time elapsed %s)\n", ts(), format(Sys.time() - t0)));
