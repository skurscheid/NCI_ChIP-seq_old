# Load libraries
library(Biostrings);
library(ggplot2);
library(ChIPsim);

# Set seed for reproducibility
set.seed(1);

# Read fasta file of reference
genome <- readDNAStringSet("GRCh38_chr22.fa.gz");

# Simulate reads assuming equal nucleotide probabilities
simRead <- function(N = 100, length = 50, alphabet = c("A", "C", "G", "T")) {
    reads <- replicate(N, paste(sample(alphabet, 50, replace = TRUE),
                                collapse = ""));
    return(DNAStringSet(reads));
}

randomQuality <- function(read) {
    paste(sample(rangePhred64, nchar(read), replace = TRUE), collapse = "");
}


# Define read quality function using Phred+64 encoding
unifQualPhred64 <- function(reads, ...) {
    rangePhred64 <- rawToChar(as.raw(64:104));
    rangePhred64 <- unlist(strsplit(rangePhred64, ""));
    qual <- lapply(reads, randomQuality);
    return(qual);
}

# Convert from specific FASTQ Phred encoding to Phred score
getReadQual <- function(qual, encoding = c("Phred33", "Phred64")) {
    score <- lapply(qual, function(x) as.numeric(charToRaw(x)));
    encoding <- match.arg(encoding);
    if (encoding == "Phred33") score <- lapply(score, function(x) x - 33);
    if (encoding == "Phred64") score <- lapply(score, function(x) x - 64);
    return(score);
}

# Plot Phred score distribution
plotReadQual <- function(qual, encoding = c("Phred33", "Phred64")) {
    score <- getReadQual(qual, encoding);
    stopifnot(length(unique(sapply(score, length))) == 1);
    df <- cbind.data.frame(pos = seq(1, unique(sapply(score, length))),
                           meanScore = rowMeans(as.data.frame(score)));
    gg <- ggplot(df, aes(x = pos, y = meanScore));
    gg <- gg + geom_line();
    gg <- gg + theme_bw();
    gg <- gg + labs(x = "Position along read",
                    y = "mean Phred score",
                    title = sprintf("N = %i reads", length(score)));
    gg;
}

reads <- simRead();
#reads[2] <- "TATAGAGCGGGTCATCGAAA";
qual <- unifQualPhred64(reads);
score <- getReadQual(qual, "Phred64");
plotReadQual(qual, "Phred64");


dfReads <- function(readPos, readNames, sequence, readLen, ...){

	## create vector to hold read sequences and qualities
	readSeq <- character(sum(sapply(readPos, sapply, length)))
	readQual <- character(sum(sapply(readPos, sapply, length)))

	idx <- 1
	## process read positions for each chromosome and strand
	for(k in length(readPos)){ ## chromosome
		for(i in 1:2){ ## strand
			for(j in 1:length(readPos[[k]][[i]])){
				## get (true) sequence
				readSeq[idx] <- as.character(readSequence(readPos[[k]][[i]][j], sequence[[k]],
						strand=ifelse(i==1, 1, -1), readLen=readLen))
				## get quality
				readQual[idx] <- randomQuality(readSeq[idx])
				## introduce sequencing errors
				readSeq[idx] <- readError(readSeq[idx], decodeQuality(readQual[idx]))
				idx <- idx + 1
			}
		}
	}
	data.frame(name=unlist(readNames), sequence=readSeq, quality=readQual,
			stringsAsFactors = FALSE)
}




myFunctions <- defaultFunctions()
myFunctions$readSequence <- dfReads
nReads <- 1000
simulated <- simChIP(nReads,
                     genome,
                     file = "output/sim_1k",
                     functions = myFunctions,
                     control = defaultControl(readDensity=list(meanLength = 150)));
