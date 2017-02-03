library(ChIPsim);
library(actuar);    # For pareto distribution

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


# Start Markov chain with background
init <- c(Binding = 0, Background = 1);
class(init) <- "StateDistribution";

# Define function to generate the parameters for background regions.
# Here we use a gamma distribution to model the background sampling weight
# for each region.
backgroundFeature <- function(
    start,
    length = 1000,
    shape = 1,
    scale = 20) {
    weight <- rgamma(1, shape = 1, scale = 20);
    params <- list(start = start, length = length, weight = weight);
    class(params) <- c("Background", "SimulatedFeature");
    return(params);
}

# Define function to generate the parameters for binding regions.
# Here we set the sampling weight for binding sites to be the average
# weight of background regions multiplied by and enrichment coefficient.
# w_binding' = t x mean(w_background). We use a Pareto distribution
# with parameter r to determine w_binding for each binding site. The
# minimum of the distribution is chosen such that its mean is w_binding'.
bindingFeature <- function(
    start,
    length = 500,
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

# Generate feature sequences
# Note that this generates a list of background and binding sites, based
# on the transition probabilities of the Markov chain model, and on the
# sampling weights for background and binding sites defined by the two
# functions backgroundFeature() and bindingFeatures().
generator <- list(
    Binding = bindingFeature,
    Background = backgroundFeature);
features <- ChIPsim::makeFeatures(
    generator,
    transition,
    init,
    start = 0,
    length = 1e6,
    globals = list(shape = 1, scale = 20),
    experimentType = "TFExperiment",
    lastFeat = c(Binding = FALSE, Background = TRUE));


# We plot the distribution of sampling weights for the derived background
# and binding sites.
df <- cbind.data.frame(
    weight = sapply(features, "[[", "weight"),
    type = sapply(features, function(x) class(x)[1]));
gg <- ggplot(df, aes(weight, fill = type));
gg <- gg + geom_density(alpha = 0.3);
gg <- gg + theme_bw();
gg <- gg + scale_fill_brewer(palette = "Dark2");
gg <- gg + xlim(0, 500);
gg <- gg + labs(
    x = "Sampling weight",
    y = "Density",
    fill = "Binding site type");
gg;

# Binding site density
constRegion <- function(weight, length) rep(weight, length);
featureDensity.Binding <- function(feature, ...) {
    constRegion(feature$weight, feature$length);
}
featureDensity.Background <- function(feature, ...) {
    constRegion(feature$weight, feature$length);
}
dens <- ChIPsim::feat2dens(features, length = 1e6);


featureDens <- lapply(features, featureDensity, background = TRUE);
