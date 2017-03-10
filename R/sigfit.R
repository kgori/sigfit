#' Fetches COSMIC's estimated mutational signatures
#' @export
fetch_cosmic_data <- function() {
    cosmic.sigs <- read.table('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt', 
                              header = TRUE, sep = '\t', check.names = FALSE)
    rownames(cosmic.sigs) <- cosmic.sigs[["Somatic Mutation Type"]]
    cosmic.sigs[, paste("Signature", 1:30)]
}

gen_bar_plot <- function(samples, featurename, title, prob, thresh) {
    feature <- rstan::extract(samples, pars = featurename)[[featurename]]
    mean_feature <- colMeans(feature)
    names(mean_feature) <- 1:length(mean_feature)
    error <- HPDinterval(as.mcmc(feature), prob = prob)
    bars <- barplot(mean_feature, ylim = c(0, max(error[, 2]*1.05)), main = title,
                    col = ifelse(error[, 1] > thresh, "dodgerblue3", "grey90"))
    top_arr <- arrows(bars, error[, 1], bars, mean_feature, angle=90, code=1, length=0.05)
    bottom_arr <- arrows(bars, error[, 2], bars, mean_feature, angle=90, code=1, length=0.05)
    list(bars, top = top_arr, bottom = bottom_arr)
}

#' Plots estimated exposure of each signature
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @importFrom "rstan" summary
#' @importFrom "rstan" extract
#' @export
plot_exposures <- function(samples, prob = 0.9, thresh = 1e-3, title = "Signature exposures") {
    plt <- gen_bar_plot(samples, "exposures", title, prob, thresh)
    plt$bars
    plt$top
    plt$bottom
    legend("topright", legend = sprintf("%.0f%% HPD > %.3f", prob*100, thresh), fill = "dodgerblue3")
}

#' Plots the fitted spectrum
#' @importFrom "rstan" extract
#' @export
plot_spectrum <- function(samples, prob = 0.9, title = "Fitted spectrum") {
    plt <- gen_bar_plot(samples, "probs", title, prob, 0)
    plt$bars
    plt$top
    plt$bottom
}

#' @useDynLib sigfit, .registration = TRUE 
#' @importFrom "rstan" sampling
#' @export
run_sampling <- function(counts, signatures, prior = NULL, ...) {
    if (is.null(prior)) {
        prior = rep(1, ncol(signatures))
    }
    stopifnot(length(counts) == nrow(signatures))
    stopifnot(length(prior) == ncol(signatures))
    dat = list(
        C = length(counts),
        S = ncol(signatures),
        counts = as.vector(counts),
        signatures = as.matrix(signatures),
        alpha = prior
    )
    model <- stanmodels$sigfit
    rstan::sampling(model, data = dat, ...)
}