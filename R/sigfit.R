#' Deal with zero values in a signature
remove_zeros_ <- function(mtx, min_allowed = 1e-9) {
    apply(mtx, 2, function(col) {
        col[col < min_allowed] <- min_allowed
        col / sum(col)
    })
}

#' Fetches COSMIC's estimated mutational signatures
#' @param reorder Reorders the matrix by substitution type and trinucleotide
#' @export
fetch_cosmic_data <- function(reorder = TRUE, remove_zeros = TRUE) {
    cosmic.sigs <- read.table('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt', 
                              header = TRUE, sep = '\t', check.names = FALSE)
    if (reorder) {
        cosmic.sigs <- cosmic.sigs[order(cosmic.sigs[["Substitution Type"]], cosmic.sigs[["Trinucleotide"]]),]
    }
    rownames(cosmic.sigs) <- cosmic.sigs[["Somatic Mutation Type"]]
    cosmic.sigs <- cosmic.sigs[, paste("Signature", 1:30)]
    if(remove_zeros) cosmic.sigs <- remove_zeros_(cosmic.sigs)
    cosmic.sigs
}

#' Access stan models
#' @export
stan_models <- function() {
    stanmodels
}

gen_bar_plot <- function(samples, featurename, title, prob, thresh, ...) {
    feature <- rstan::extract(samples, pars = featurename)[[featurename]]
    mean_feature <- colMeans(feature)
    names(mean_feature) <- 1:length(mean_feature)
    error <- HPDinterval(as.mcmc(feature), prob = prob)
    bars <- barplot(mean_feature, ylim = c(0, max(error[, 2]*1.05)), main = title,
                    col = ifelse(error[, 1] > thresh, "dodgerblue3", "grey90"), ...)
    top_arr <- arrows(bars, error[, 1], bars, mean_feature, angle=90, code=1, length=0.05)
    bottom_arr <- arrows(bars, error[, 2], bars, mean_feature, angle=90, code=1, length=0.05)
    list(bars, top = top_arr, bottom = bottom_arr)
}

#' Plots estimated exposure of each signature
#' @param samples The MCMC samples
#' @param prob The width of the HPD interval
#' @param thresh Signatures with a lower HPDI below this value are coloured grey
#' @param title The main title of this plot
#' @param ... Arguments to pass through to graphics::barplot
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @importFrom "rstan" summary
#' @importFrom "rstan" extract
#' @export
plot_exposures <- function(samples, prob = 0.9, thresh = 1e-3, title = "Signature exposures", ...) {
    plt <- gen_bar_plot(samples, "exposures", title, prob, thresh, ...)
    plt$bars
    plt$top
    plt$bottom
    legend("topright", legend = sprintf("%.0f%% HPD > %.3f", prob*100, thresh), fill = "dodgerblue3")
}

#' Plots the fitted spectrum
#' @param samples The MCMC samples
#' @param prob The width of the HPD interval
#' @param title The main title of this plot
#' @importFrom "rstan" extract
#' @export
plot_spectrum <- function(samples, prob = 0.9, title = "Fitted spectrum", ...) {
    plt <- gen_bar_plot(samples, "probs", title, prob, 0, ...)
    plt$bars
    plt$top
    plt$bottom
}

#' Runs the MCMC sampling chain to estimate exposures
#' @param counts Vector of mutation counts
#' @param signatures Matrix of mutational signatures
#' @param prior Vector of the same length as signatures, to be used as the Dirichlet prior in the sampling chain. Default is all ones, i.e. uninformative
#' @param ... Arguments to pass to rstan::sampling
#' @examples
#'  # Custom prior favours signature 1 over 2, 3 and 4
#' samples <- sigfit::run_sampling(mycounts, mysignatures, prior = c(5, 1, 1, 1))
#' 
#' # Run a single chain for quite a long time
#' samples <- sigfit::run_sampling(mycounts, mysignatures, chains = 1, niter = 13000, warmup = 3000)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @export
run_sampling <- function(counts, signatures, prior = NULL, hierarchical = FALSE, multi = FALSE, ...) {
    if (is.null(prior)) {
        prior = rep(1, ncol(signatures))
    }
    
    # If counts is a vector, convert it to a single row matrix
    if (is.vector(counts)) counts <- matrix(counts, nrow = 1)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NCAT, NSIG]
    stopifnot(ncol(counts) == nrow(signatures))
    stopifnot(length(prior) == ncol(signatures))
    
    dat = list(
        C = ncol(counts),
        S = ncol(signatures),
        G = nrow(counts),
        counts = as.matrix(counts),
        signatures = as.matrix(signatures),
        alpha = prior
    )
    if (hierarchical) {
        model <- stanmodels$sigfit_hier
    }
    else {
        if (multi) {
            model <- stanmodels$sigfit_multi
        }
        else {
            model <- stanmodels$sigfit
        }
    }
    rstan::sampling(model, data = dat, ...)
}

#' Extracts signatures from a set of mutation counts
#' using models based on NMF or EMu
#' 
#' @param counts Matrix of mutation counts for each sample (rows) in each category
#' (columns)
#' @param nsignatures Number of signatures to extract
#' @param method Either "emu" or "nmf" (though currently "nmf" is experimental)
#' @param opportunities Optional matrix of mutational opportunities for "emu" method.
#' Dimensions should be same as for counts
#' @param stanfunc "sampling"|"optimizing"|"vb" Choice of rstan inference strategy. 
#' "sampling" is the full Bayesian MCMC approach, and is the default. "optimizing"
#' returns the Maximum a Posteriori (MAP) point estimates via numerical optimization.
#' "vb" uses Variational Bayes to approximate the full posterior.
#' @param ... Any other parameters to pass through to rstan 
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @export
extract_signatures <- function(counts, nsignatures, method = "emu", 
                               opportunities = NULL, exposures_prior = 0.5, 
                               stanfunc = "sampling", ...) {
    if (method == "emu") {
        if(is.null(opportunities)) {
            opportunities <- matrix(1, nrow = nrow(counts), ncol = ncol(counts))
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        model <- stanmodels$sigfit_emu
        data <- list(
            N = ncol(counts),
            M = nrow(counts),
            n = nsignatures,
            counts = as.matrix(counts),
            opps = as.matrix(opportunities)
        )
    }
    
    else if (method == "nmf") {
        model <- stanmodels$sigfit_nmf
        data <- list(
            G = nrow(counts),
            C = ncol(counts),
            S = nsignatures,
            counts = counts,
            exposures_prior_val = exposures_prior
        )
    }
    
    if (stanfunc == "sampling") {
        cat("Stan sampling:")
        return(sampling(model, data = data, chains = 1, ...))
    }
    else if (stanfunc == "optimizing") {
        cat("Stan optimizing:")
        return(optimizing(model, data = data, ...))
    }
    else if (stanfunc == "vb") {
        cat("Stan vb")
        return(vb(model, data = data, ...))
    }
}
