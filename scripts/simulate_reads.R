#!/usr/bin/env Rscript

# The R command line script simulates reads for a TF ChIP experiment using
# the R package ChIPsim by Peter Humburg
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 14/02/2017
# Last change: 16/02/2017


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
        c("--simName"),
        type = "character",
        default = "test",
        help = "Identifier of simulation [default %default]",
        metavar = "character"),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        default = "rawData",
        help = "Output directory [default %default]",
        metavar = "character"),
    make_option(
        c("--bindProb"),
        type = "double",
        default = 0.05,
        help = "P(binding|background) in MC model [default %default]",
        metavar = "double"),
    make_option(
        c("--backProb"),
        type = "double",
        default = 0.95,
        help = "P(background|background) in MC model [default %default]",
        metavar = "double"),
    make_option(
        c("--EF"),
        type = "double",
        default = 5,
        help = "Enrichment factor of mean(TF binding)/mean(Background binding) [default %default]",
        metavar = "double"),
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
        metavar = "integer"),
    make_option(
        c("--gSize"),
        type = "integer",
        default = 0,
        help = "Genome size [in bp]; if not =0, then the number of reads is scaled such that nReads = nReads/gSize * size(ref) [default %default]",
        metavar = "integer"),
    make_option(
        c("--nReps"),
        type = "integer",
        default = 3,
        help = "Number of replicates [default %default]",
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
simName <- args$simName;
outdir <- args$outdir;
Pbind_given_back <- args$bindProb;
Pback_given_back <- args$backProb;
EF <- args$EF;
seed <- args$seed;
backgroundLength <- args$backLength;
bindingLength <- args$bindLength;
meanFragmentLength <- args$meanFragLength;
minFragmentLength <- args$minFragLength;
maxFragmentLength <- args$maxFragLength;
readLength <- args$readLength;
nReads <- args$nReads;
gSize <- args$gSize;
nReps <- args$nReps;
cat(sprintf("%s Parameter summary\n", ts()));
cat(sprintf(" Reference          = %s\n", genome));
cat(sprintf(" name               = %s\n", simName));
cat(sprintf(" outdir             = %s\n", outdir));
cat(sprintf(" seed               = %s\n", seed));
cat(sprintf(" Pbind_given_back   = %f\n", Pbind_given_back));
cat(sprintf(" Pback_given_back   = %f\n", Pback_given_back));
cat(sprintf(" EF                 = %f\n", EF));
cat(sprintf(" bindingLength      = %i\n", bindingLength));
cat(sprintf(" backgroundLength   = %i\n", backgroundLength));
cat(sprintf(" meanFragmentLength = %i\n", meanFragmentLength));
cat(sprintf(" minFragmentLength  = %i\n", minFragmentLength));
cat(sprintf(" maxFragmentLength  = %i\n", maxFragmentLength));
cat(sprintf(" readLength         = %i\n", readLength));
cat(sprintf(" nReads             = %i\n", nReads));
cat(sprintf(" gSize              = %i\n", gSize));
cat(sprintf(" nReps              = %i\n", nReps));


## ------------------------------------------------------------------------
## Read in reference
genome <- readDNAStringSet(genome);
refSize <- width(genome);
cat(sprintf(
    "%s Read reference sequence \"%s\" of length %i\n",
    ts(),
    names(genome),
    refSize));
if (refSize < backgroundLength + bindingLength) {
    stop("Reference sequence must be longer than bindingLength + backgroundLength.\n", call. = FALSE);
}
if (gSize > 0 & gSize < refSize) {
    stop(sprintf("gSize = %i < refSize = %i!", gSize, refSize));
}



## ------------------------------------------------------------------------
## Define transition probabilities and initial state
transition <- list(
    Binding = c(Background = 1),
    Background = c(Binding = Pbind_given_back, Background = Pback_given_back)
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
    length = refSize,
    globals = list(shape = 1, scale = 20),
    control = list(Binding = list(enrichment = EF)),
    experimentType = "TFExperiment",
    lastFeat = c(Binding = FALSE, Background = TRUE));


## ------------------------------------------------------------------------
## Stop if less <=2 binding sites
t <- table(sapply(features, class));
nBinding <- t[which(names(t) == "Binding")];
nBinding <- ifelse(length(nBinding) == 0, 0, nBinding);
nBackground <- t[which(names(t) == "Background")];
nBackground <- ifelse(length(nBackground) == 0, 0, nBackground);
cat(sprintf("%s Results of MC model\n", ts()));
cat(sprintf(" %10i binding sites\n", nBinding));
cat(sprintf(" %10i background sites\n", nBackground));
if (nBinding < 1 & Pbind_given_back > 0) {
    stop("MC model gave <1 binding sites. Decrease binding site lengths or change reference.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Store binding sites in BED file
cat(sprintf("%s Storing TF binding sites in BED file...\n", ts()));
chr <- unlist(strsplit(names(genome), " "))[1];
bindFeat <- features[sapply(features, class)[1, ] == "Binding"];
df.BED <- cbind.data.frame(
    chr = rep(chr, length(bindFeat)),
    start = sapply(bindFeat, "[[", 1) - 1,
    end = sapply(bindFeat, "[[", 1) + sapply(bindFeat, "[[", 2),
    name = sprintf("bindFeat_%s_%i", simName, seq(1, length(bindFeat))),
    score = sapply(bindFeat, "[[", 3),
    strand = rep(".", length(bindFeat)));
write.table(
    df.BED,
    file = sprintf("%s/bindingSites_%s.bed", outdir, simName),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE);


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
    filename = sprintf(
        "%s/distr_binding_site_strength_%s.pdf",
        outdir,
        simName),
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
dens <- ChIPsim::feat2dens(features, length = refSize);


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
## Define uniform read quality function
randomQualityPhred33 <- function(read, ...) {
    # Character vector of symbols for the Phred+33 quality encoding scale
    rangePhred33 <- unlist(strsplit(rawToChar(as.raw(33:126)), ""));
    # Uniform-randomly sample qualities
    paste(sample(rangePhred33, width(read), replace = TRUE), collapse = "");
}


## ------------------------------------------------------------------------
## Sample reads for every replicate
if (gSize != 0) {
    nReads <- round(nReads / gSize * refSize);
}
cat(sprintf("%s Sampling %i reads for %i replicates...\n", ts(), nReads, nReps));
for (i in 1:nReps) {
    # Sample read/fragment starting position
    fragStart <- ChIPsim::sampleReads(readDens, nreads = nReads);

    # Sample fragment lengths and store start and length in list
    frag <- lapply(fragStart, function(x) {
        cbind.data.frame(
            start = x,
            length = sample(
                minFragmentLength:maxFragmentLength,
                length(x),
                replace = TRUE,
                prob = fragLength(
                    minFragmentLength:maxFragmentLength,
                    minLength = minFragmentLength,
                    maxLength = maxFragmentLength,
                    meanLength = meanFragmentLength)))
    });

    # We need to make sure that
    #    (1) fragstart + fraglength <= refSize on the positive strand
    #    (2) fragstart - fraglength >= 1 on the negative strand
    idx1 <- frag[[1]]$start + frag[[1]]$length - 1 <= refSize;
    idx2 <- frag[[2]]$start - frag[[2]]$length + 1 >= 1;
    frag[[1]] <- frag[[1]][idx1, ];
    frag[[2]] <- frag[[2]][idx2, ];

    # Extract sequences and write to two FASTQ
    cat(sprintf("%s Writing pe reads of rep %i to fastq files...\n", ts(), i));
    # Open two files for the left and right reads
    file.R1 <- gzfile(
        sprintf(
            "%s/reads_%s_rep%i_R1.fastq.gz",
            outdir,
            simName,
            i),
        "w");
    file.R2 <- gzfile(
        sprintf(
            "%s/reads_%s_rep%i_R2.fastq.gz",
            outdir,
            simName,
            i),
        "w");
    for (j in 1:2) {
        for (k in 1:nrow(frag[[j]])) {
            if (j == 1) {
                # frag[j==1] gives 5' start positions from the positive strand
                #                                read2
                #                             <xxxxxx|
                # rev = 3' oooooooooooooooooooooooooooooooooooooooooooooo 5'
                # frag        |----------------------|
                # fragstart   x
                # fwd = 5' oooooooooooooooooooooooooooooooooooooooooooooo 3'
                #             |xxxxxx>
                #             read1
                read1 <- subseq(
                    genome,
                    start = frag[[j]]$start[k],
                    end = frag[[j]]$start[k] + readLength - 1);
                read2 <- subseq(
                    genome,
                    start = frag[[j]]$start[k] + frag[[j]]$length[k] - readLength,
                    end = frag[[j]]$start[k] + frag[[j]]$length[k] - 1);
                read2 <- reverseComplement(read2);
            } else {
                # frag[j==2] gives 5' start positions from the negative strand
                #                                read1
                #                             <xxxxxx|
                # rev = 3' oooooooooooooooooooooooooooooooooooooooooooooo 5'
                # frag        |----------------------|
                # fragstart                          x
                # fwd = 5' oooooooooooooooooooooooooooooooooooooooooooooo 3'
                #             |xxxxxx>
                #             read2
                read1 <- subseq(
                    genome,
                    start = frag[[j]]$start[k] - readLength + 1,
                    end = frag[[j]]$start[k]);
                read1 <- reverseComplement(read1);
                read2 <- subseq(
                    genome,
                    start = frag[[j]]$start[k] - frag[[j]]$length[k] + 1,
                    end = frag[[j]]$start[k] - frag[[j]]$length[k] + readLength);
            }
            id <- sprintf("simul_%s_rep%i_strand%i_read%i", simName, i, j, k);
            # Write left read
            cat("@", id, "/1\n", sep = "", file = file.R1);
            cat(as.character(read1), "\n", sep ="", file = file.R1);
            cat("+", "\n", sep = "", file = file.R1);
            cat(randomQualityPhred33(read1), "\n", sep = "", file = file.R1);
            # Write right read
            cat("@", id, "/2\n", sep = "", file = file.R2);
            cat(as.character(read2), "\n", sep ="", file = file.R2);
            cat("+", "\n", sep = "", file = file.R2);
            cat(randomQualityPhred33(read2), "\n", sep = "", file = file.R2);
        }
    }
    close(file.R1);
    close(file.R2);
}


## ------------------------------------------------------------------------
## Done
runTime <- paste("time elapsed", format(Sys.time() - t0));
cat(sprintf("%s Done (%s)\n", ts(), runTime));
