## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.path = 'figures/', echo = TRUE)

## ------------------------------------------------------------------------
suppressMessages(library(Biostrings));  # For DNAStringSet
suppressMessages(library(ChIPsim));     # Main package for ChIP-seq data simulation
suppressMessages(library(ggplot2));     # For plotting
suppressMessages(library(actuar));      # For pareto distribution
suppressMessages(library(reshape2));    # For melt (convert wide df to long df)
set.seed(1234);

## ------------------------------------------------------------------------
genome <- readDNAStringSet("GRCh38_MT.fa.gz");
refLength <- width(genome);
print(refLength);

## ------------------------------------------------------------------------
# Transition probabilities
# Note that a binding region has to be followed by a background region,
# while background regions are allowed to follow each other.
#       Bind  Back
#       ----------
# P =   0.00  1.00  | Bind
#       0.05  0.95  | Back
transition <- list(
    Binding = c(Background = 1),
    Background = c(Binding = 0.05, Background = 0.95)
);
transition <- lapply(transition, "class<-", "StateDistribution");

## ------------------------------------------------------------------------
# Start Markov chain with background
init <- c(Binding = 0, Background = 1);
class(init) <- "StateDistribution";

## ------------------------------------------------------------------------
# Define background region length
#backgroundLength <- 1000;
backgroundLength <- 500;
# Define function to generate the parameters for background regions.
# Here we use a gamma distribution to model the background sampling weight
# for each region.
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
# Define binding region length
bindingLength = 50;
# Define function to generate the parameters for binding regions.
# Here we set the sampling weight for binding sites to be the average
# weight of background regions multiplied by and enrichment coefficient.
# w_binding' = t x mean(w_background). We use a Pareto distribution
# with parameter r to determine w_binding for each binding site. The
# minimum of the distribution is chosen such that its mean is w_binding'.
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
# Generate feature sequences
# Note that this generates a list of background and binding sites, based
# on the transition probabilities of the Markov chain model, and on the
# sampling weights for background and binding sites defined by the two
# functions backgroundFeature() and bindingFeatures().
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
table(sapply(features, class));

## ------------------------------------------------------------------------
features[[1]];

## ------------------------------------------------------------------------
# We plot the distribution of sampling weights for the derived background
# and binding sites.
df <- cbind.data.frame(
    weight = sapply(features, "[[", "weight"),
    type = sapply(features, function(x) class(x)[1]));
gg <- ggplot(subset(df, weight <= 500), aes(weight, fill = type));
gg <- gg + geom_density(alpha = 0.3);
gg <- gg + theme_bw();
gg <- gg + scale_fill_brewer(palette = "Dark2");
gg <- gg + xlim(0, 500);
gg <- gg + labs(
    x = "Sampling weight",
    y = "Density",
    fill = "Binding site type");
gg;

## ------------------------------------------------------------------------
# Binding site density
featureDensity.Binding <- function(feature, ...) {
  featDens <- numeric(feature$length);
  featDens[floor(feature$length / 2)] <- feature$weight;
#  featDens[floor(feature$length / 2)] <- feature$weight * feature$length;
  return(featDens);
}

## ------------------------------------------------------------------------
featureDensity.Background <- function(feature, ...) {
  featDens <- numeric(feature$length);
#  featDens[] <- feature$weight;
  featDens[] <- feature$weight / bindingLength;
  return(featDens);
}

## ------------------------------------------------------------------------
dens <- ChIPsim::feat2dens(features, length = refLength);

## ------------------------------------------------------------------------
df.site <- cbind.data.frame(
  pos = seq(1, length(dens)), 
  dens = dens);
df2 <- cbind.data.frame(
  pos = sapply(features, "[[", 1), 
  weight = sapply(features, "[[", 3), 
  sapply(features, class)[1, ]);
gg <- ggplot(df.site, aes(x = pos, y = dens));
gg <- gg + geom_line();
gg <- gg + geom_step(data = df2, aes(x = pos, y = weight), colour = "red", alpha = 0.5);
gg <- gg + theme_bw();
gg <- gg + labs(
  x = "Position along reference",
  y = "Density");
gg;

## ------------------------------------------------------------------------
range <- 5000;
xpeak <- which.max(df.site$dens);   # Position of largest peak
gg <- ggplot(df.site[xpeak - seq(-range, range), ]);
gg <- gg + geom_line(aes(x = pos, y = dens));
gg <- gg + geom_step(data = df2[which(abs(df2$pos - xpeak) <= range), ], 
                     aes(x = pos, y = weight), 
                     colour = "red", 
                     alpha = 0.5);
gg <- gg + theme_bw();
gg <- gg + labs(
  x = "Position along reference",
  y = "Density");
gg <- gg + scale_y_log10();
gg;

## ------------------------------------------------------------------------
fragLength <- function(x, minLength, maxLength, meanLength, ...) {
  sd <- (maxLength - minLength)/4;
  prob <- dnorm(minLength:maxLength, mean = meanLength, sd = sd);
  prob <- prob/sum(prob);
  return(prob[x - minLength + 1]);
}

## ------------------------------------------------------------------------
readDens <- ChIPsim::bindDens2readDens(
  dens, 
  fragLength, 
  bind = 50, 
  minLength = 150, 
  maxLength = 250,
  meanLength = 200);

## ------------------------------------------------------------------------
range <- 500;
df.DNA <- cbind.data.frame(
  pos = seq(1, nrow(readDens)),
  readDens);
colnames(df.DNA)[2:3] <- c("positive", "negative");
df.DNA <- melt(df.DNA, id.vars = "pos");
df <- rbind.data.frame(
  cbind.data.frame(pos = df.site$pos, 
                   variable = "bindingSite", 
                   value = df.site$dens),
  df.DNA);
gg <- ggplot(subset(df, abs(pos - xpeak) <= range));
gg <- gg + geom_line(aes(x = pos, y = value));
gg <- gg + facet_grid(variable ~ ., scales = "free_y");
gg <- gg + theme_bw();
gg <- gg + labs(
  x = "Position along reference",
  y = "Density");
gg <- gg + scale_y_log10();
gg;

## ------------------------------------------------------------------------
readLoc <- ChIPsim::sampleReads(readDens, nreads = 1e5);

## ------------------------------------------------------------------------
randomQualityPhred33 <- function(read, ...) {
  # Character vector of symbols for the Phred+33 quality encoding scale
  rangePhred33 <- unlist(strsplit(rawToChar(as.raw(33:126)), ""));
  # Uniform-randomly sample qualities
  paste(sample(rangePhred33, length(read), replace = TRUE), collapse = "");
}

## ------------------------------------------------------------------------
# Read length
readLength <- 100;

# We need to make sure that readLoc + readLen <= refLength for both strands
readLoc[[1]] <- readLoc[[1]][which(readLoc[[1]] + readLength <= refLength)];
readLoc[[2]] <- readLoc[[1]][which(readLoc[[1]] - readLength > 0)];

# Create names
nreads <- sapply(readLoc, length);
names <- list(fwd = sprintf("read_fwd_%s", seq(nreads[1])),
              rev = sprintf("read_rev_%s", seq(nreads[2])));

# Write to FASTQ
# Uncomment for output
#pos2fastq(readLoc, 
#          names = names,
#          sequence = genome[[1]],
#          qualityFun = randomQualityPhred33, 
#          errorFun = readError,
#          readLen = readLength,
#          file = "TF_ChIP_MT.fastq");

