#!/usr/bin/env Rscript

# The R command line script reads in fa/fa.gz files from a directory
# and stores the size of each file in a summary text files
# This is basically an R implementation of "samtools faidx" using the R
# package Biostrings; "samtools faidx" unfortunately does not allow indexing
# gzipped fasta files
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 21/02/2017
# Last change: 21/02/2017


## ------------------------------------------------------------------------
## Start clock
t0 <- Sys.time();


## ------------------------------------------------------------------------
## Load libraries
suppressMessages(library(optparse));    # For Python-style command line args
suppressMessages(library(Biostrings));  # For DNAStringSet


## ------------------------------------------------------------------------
## Parse command line arguments
option_list <- list(
    make_option(
        c("-d", "--dir"),
        type = "character",
        default = NULL,
        help = "The directory containing the fa/fa.gz files",
        metavar = "character"),
    make_option(
        c("-r", "--regexp"),
        type = "character",
        default = "",
        help = "Regular expression used for filtering files [default %default]",
        metavar = "character"),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "fastq_index.csv",
        help = "Output summary text file [default %default]",
        metavar = "character")
    );
opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);
if (is.null(args$dir)) {
    print_help(opt_parser);
    stop("Please specify a directory.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Timestamp function
ts <- function() {
    return(format(Sys.time(), "[%a %b %d %Y %H:%M:%S]"));
}


## ------------------------------------------------------------------------
## Global variables
dir <- args$dir;
regexp <- args$regexp;
out <- args$out;
cat(sprintf("%s Parameter summary\n", ts()));
cat(sprintf(" dir    = %s\n", dir));
cat(sprintf(" regexp = %s\n", regexp));
cat(sprintf(" out    = %s\n", out));


## ------------------------------------------------------------------------
## List files in directory
files <- list.files(
    path = dir,
    pattern = "(fa|fa.gz|fasta|fasta.gz)",
    full.names = TRUE,
    no.. = TRUE,
    include.dirs = FALSE);
cat(sprintf("%s Found %s files in %s:\n", ts(), length(files), dir));
if (length(files) > 0) {
    for (i in 1:length(files)) {
        cat(sprintf(" %s\n", files[i]))
    }
}


## ------------------------------------------------------------------------
## Read files
DNA <- lapply(files, readDNAStringSet);


## ------------------------------------------------------------------------
## Write to output file
df <- cbind.data.frame(
    file = basename(files),
    length = sapply(DNA, width),
    name = sapply(DNA, names));
write.table(df, file = out, sep = ",", row.names = FALSE, col.names = FALSE);


## ------------------------------------------------------------------------
## Done
runTime <- paste("time elapsed", format(Sys.time() - t0));
cat(sprintf("%s Done (%s)\n", ts(), runTime));
