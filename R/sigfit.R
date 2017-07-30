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

#' Returns human genome or exome trinucleotide frequencies. This is useful to
#' de-normalize signatures that were obtained using the mutational opportunities
#' from the human genome/exome, in order to compare them to signatures extracted
#' without incorporating opportunities.
#' @param type Either "genome" (default) or "exome".
#' @examples
#' Extract signatures using human exome opportunitites
#' samples <- sigfit::extract_signatures(mycounts, nsignatures = 3, method = "emu", opportunities = "human-exome")
#' sigs <- retrieve_pars(samples, "signatures")
#' 
#' # De-normalize (mean) extracted signatures
#' freqs <- human_trinuc_freqs("exome")
#' denorm_sigs <- apply(sigs, 1, function(sig) { 
#'      tmp <- sig[,1] * freqs
#'      tmp / sum(tmp)
#' })
#' @useDynLib sigfit, .registration = TRUE
#' @export
human_trinuc_freqs <- function(type = "genome") {
    if (type == "genome") {
        # Human genome trinucleotide frequencies (from EMu)
        matrix(rep(c(1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>A @ AC[ACGT]
                     1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, # C>A @ CC[ACGT]
                     8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, # C>A @ GC[ACGT]
                     1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08, # C>A @ TC[ACGT]
                     1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>G @ AC[ACGT]
                     1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, # C>G @ CC[ACGT]
                     8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, # C>G @ GC[ACGT]
                     1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08, # C>G @ TC[ACGT]
                     1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>T @ AC[ACGT]
                     1.05e+08, 7.46e+07, 1.57e+07, 1.01e+08, # C>T @ CC[ACGT]
                     8.17e+07, 6.76e+07, 1.35e+07, 7.93e+07, # C>T @ GC[ACGT]
                     1.11e+08, 8.75e+07, 1.25e+07, 1.25e+08, # C>T @ TC[ACGT]
                     1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, # T>A @ AC[ACGT]
                     7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, # T>A @ CC[ACGT]
                     6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07, # T>A @ GC[ACGT]
                     1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, # T>A @ TC[ACGT]
                     1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, # T>C @ AC[ACGT]
                     7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, # T>C @ CC[ACGT]
                     6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07, # T>C @ GC[ACGT]
                     1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, # T>C @ TC[ACGT]
                     1.17e+08, 7.57e+07, 1.04e+08, 1.41e+08, # T>G @ AC[ACGT]
                     7.31e+07, 9.55e+07, 1.15e+08, 1.13e+08, # T>G @ AC[ACGT]
                     6.43e+07, 5.36e+07, 8.52e+07, 8.27e+07, # T>G @ AG[ACGT]
                     1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08),# T>G @ AT[ACGT]
                   nrow(counts)),
               nrow = nrow(counts), ncol = ncol(counts), byrow = T)
    }
    else if (type == "exome") {
        # Human exome trinucleotide frequencies (from EMu)
        matrix(rep(c(1940794, 1442408, 514826, 1403756,
                     2277398, 2318284, 774498, 2269674,
                     1740752, 1968596, 631872, 1734468,
                     1799540, 1910984, 398440, 2024770,
                     1940794, 1442408, 514826, 1403756,
                     2277398, 2318284, 774498, 2269674,
                     1740752, 1968596, 631872, 1734468,
                     1799540, 1910984, 398440, 2024770,
                     1940794, 1442408, 514826, 1403756,
                     2277398, 2318284, 774498, 2269674,
                     1740752, 1968596, 631872, 1734468,
                     1799540, 1910984, 398440, 2024770,
                     1299256, 1166912, 1555012, 1689928,
                     978400,  2119248, 2650754, 1684488,
                     884052,  1173252, 1993110, 1251508,
                     1391660, 1674368, 1559846, 2850934,
                     1299256, 1166912, 1555012, 1689928,
                     978400,  2119248, 2650754, 1684488,
                     884052,  1173252, 1993110, 1251508,
                     1391660, 1674368, 1559846, 2850934,
                     1299256, 1166912, 1555012, 1689928,
                     978400,  2119248, 2650754, 1684488,
                     884052,  1173252, 1993110, 1251508,
                     1391660, 1674368, 1559846, 2850934),
                   nrow(counts)),
               nrow = nrow(counts), ncol = ncol(counts), byrow = T)
    }
    else {
        stop("type must be either \"genome\" or \"exome\"")
    }
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

#' Obtains summary values for a set of model parameters (signatures or exposures) from a stanfit object.
#' @param object An object of class stanfit.
#' @param feature Name of the parameter set to extract; either "signatures" or "exposures".
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the HPD intervals.
#' @param signature_names Vector containing the names of the signatures used for fitting. Used only when 
#' retrieving exposures from fitted signatures.
#' @examples
#' # Extract signatures using the EMu (Poisson) model
#' samples <- sigfit::extract_signatures(mycounts, nsignatures = 3, method = "emu", opportunities = "human-genome")
#' 
#' # Retrieve array of signatures
#' signatures <- sigfit::retrieve_pars(samples, "signatures")
#' 
#' # Retrieve array of exposures using 90% HPD intervals
#' exposures <- sigfit::retrieve_pars(samples, "exposures", prob = 0.9)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
#' @export
retrieve_pars <- function(object, feature, prob = 0.95, signature_names = NULL) {
    feat <- rstan::extract(object, pars = feature)[[feature]]
    # Multi-sample case
    if (length(dim(feat)) > 2) {
        names1 <- ifelse(feature == "signatures",
                         ifelse(!is.null(signature_names),
                                signature_names,
                                paste("Signature", LETTERS[1:dim(feat)[2]])),
                         NULL)
        names2 <- ifelse(feature == "exposures",
                         ifelse(!is.null(signature_names),
                                signature_names,
                                paste("Signature", LETTERS[1:dim(feat)[3]])),
                         NULL)
        # for signatures: dims = (signatures, categories, mean/lwr/upr)
        # for exposures: dims = (samples, signatures, mean/lwr/upr)
        feat.summ <- array(NA, dim = c(dim(feat)[2], dim(feat)[3], 3),
                           dimnames = list(names1, names2, c("mean", paste0(c("lower_", "upper_"), prob))))
        for (i in 1:dim(feat.summ)[1]) {
            feat.summ[i,,1] <- colMeans(feat[,i,])
            feat.summ[i,,2:3] <- coda::HPDinterval(coda::as.mcmc(feat[,i,]), prob = prob)
        }
    } 
    # Single-sample case
    else {
        names1 <- ifelse(feature == "signatures",
                         ifelse(!is.null(signature_names),
                                signature_names,
                                "Signature A"),
                         NULL)
        names2 <- ifelse(feature == "exposures",
                         ifelse(!is.null(signature_names),
                                signature_names,
                                paste("Signature", LETTERS[1:dim(feat)[2]])),
                         NULL)
        feat.summ <- array(NA, dim = c(1, dim(feat)[2], 3),
                           dimnames = list(names1, names2, c("mean", paste0(c("lower_", "upper_"), prob))))
        feat.summ[1,,1] <- colMeans(feat)
        feat.summ[1,,2:3] <- coda::HPDinterval(coda::as.mcmc(feat), prob = prob)
    }
    feat.summ
}

#' Runs MCMC to fit signatures and estimate exposures
#' @param counts Matrix of mutation counts per category (columns) per genome sample (rows).
#' @param signatures Matrix of mutational signatures (columns) to be fitted.
#' @param prior Vector of the same length as signatures, to be used as the Dirichlet prior in the sampling chain. Default prior is uniform (uninformative).
#' @param method Either "emu" or "nmf".
#' @param ... Arguments to pass to rstan::sampling.
#' @examples
#'  # Custom prior favours signature 1 over 2, 3 and 4
#' samples <- sigfit::fit_signatures(mycounts, mysignatures, prior = c(5, 1, 1, 1))
#' 
#' # Run a single chain for quite a long time
#' samples <- sigfit::fit_signatures(mycounts, mysignatures, chains = 1, niter = 13000, warmup = 3000)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @export
fit_signatures <- function(counts, signatures, prior = NULL, hierarchical = FALSE, method = "nmf", ...) {
    if (is.null(prior)) {
        prior = rep(1, ncol(signatures))
    }
    
    # If counts is a vector, convert it to a single row matrix
    if (is.vector(counts)) counts <- matrix(counts, nrow = 1)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NCAT, NSIG]
    stopifnot(ncol(counts) == nrow(signatures))
    stopifnot(length(prior) == ncol(signatures))
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    if (method == "emu") {
        # NEED TO IMPLEMENT alpha
        dat = list(
            N = ncol(counts),
            n = ncol(signatures),
            M = nrow(counts),
            mu = t(as.matrix(signatures)),
            X = t(as.matrix(counts)),
            w = t(as.matrix(opportunities))
        )
        model <- stanmodels$sigfit_fit_emu
    }
    else if (hierarchical) {
        dat = list(
            C = ncol(counts),
            S = ncol(signatures),
            G = nrow(counts),
            signatures = as.matrix(signatures),
            counts = as.matrix(counts)
        )
        model <- stanmodels$sigfit_fit_nmf_hier
    }
    else {
        dat = list(
            C = ncol(counts),
            S = ncol(signatures),
            G = nrow(counts),
            signatures = as.matrix(signatures),
            counts = as.matrix(counts),
            alpha = prior
        )
        model <- stanmodels$sigfit_fit_nmf
    }
    rstan::sampling(model, data = dat, ...)
}

#' Extracts signatures from a set of mutation counts
#' using models based on NMF or EMu
#' 
#' @param counts Matrix of mutation counts for each sample (rows) in each category
#' (columns).
#' @param nsignatures Number of signatures to extract.
#' @param method Either "emu" or "nmf" (though currently "nmf" is experimental).
#' @param opportunities Optional matrix of mutational opportunities for "emu" method; must have same dimension as counts. 
#' If equals to "human-genome" or "human-exome", the reference human genome/exome opportunities will be used for every sample.
#' @param stanfunc "sampling"|"optimizing"|"vb" Choice of rstan inference strategy. 
#' "sampling" is the full Bayesian MCMC approach, and is the default. "optimizing"
#' returns the Maximum a Posteriori (MAP) point estimates via numerical optimization.
#' "vb" uses Variational Bayes to approximate the full posterior.
#' @param ... Any other parameters to pass through to rstan.
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @export
extract_signatures <- function(counts, nsignatures, method = "emu", 
                               opportunities = NULL, exposures_prior = 0.5, 
                               stanfunc = "sampling", ...) {
    if (method == "emu") {
        if (is.null(opportunities)) {
            opportunities <- matrix(1, nrow = nrow(counts), ncol = ncol(counts))
        }
        else if (opportunities == "human-genome") {
            opportunities <- human_trinuc_freqs("genome")
        }
        else if (opportunities == "human-exome") {
            opportunities <- human_trinuc_freqs("exome")
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        model <- stanmodels$sigfit_ext_emu
        data <- list(
            N = ncol(counts),
            M = nrow(counts),
            n = nsignatures,
            counts = as.matrix(counts),
            opps = as.matrix(opportunities)
        )
    }
    
    else if (method == "nmf") {
        if (!is.null(opportunities)) {
            warning("Using \"nmf\" model; opportunities will be ignored.")
        }
        model <- stanmodels$sigfit_ext_nmf
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

#' Fits signatures to estimate exposures in a set of mutation counts
#' and extracts additional signatures present in the samples.
#' 
#' @param counts Matrix of mutation counts for each sample (rows) in each category
#' (columns).
#' @param signatures Matrix of fixed mutational signatures (columns) to be fitted.
#' @param num_extra_sigs Number of additional signatures to be extracted.
#' @param stanfunc "sampling"|"optimizing"|"vb" Choice of rstan inference strategy. 
#' "sampling" is the full Bayesian MCMC approach, and is the default. "optimizing"
#' returns the Maximum a Posteriori (MAP) point estimates via numerical optimization.
#' "vb" uses Variational Bayes to approximate the full posterior.
#' @param ... Any other parameters to pass through to rstan.
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @export
fit_extract_signatures <- function(counts, signatures, num_extra_sigs, 
                                   stanfunc = "sampling", ...) {
    # If counts is a vector, convert it to a single row matrix
    if (is.vector(counts)) counts <- matrix(counts, nrow = 1)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NCAT, NSIG]
    stopifnot(ncol(counts) == nrow(signatures))
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    dat = list(
        C = ncol(counts),
        S = ncol(signatures),
        G = nrow(counts),
        N = num_extra_sigs,
        fixed_sigs = t(as.matrix(signatures)),
        counts = as.matrix(counts)
    )
    ## Only NMF implemented so far
    model <- stanmodels$sigfit_fitex_nmf
    
    if (stanfunc == "sampling") {
        cat("Stan sampling:")
        return(rstan::sampling(model, data = data, chains = 1, ...))
    }
    else if (stanfunc == "optimizing") {
        cat("Stan optimizing:")
        return(rstan::optimizing(model, data = data, ...))
    }
    else if (stanfunc == "vb") {
        cat("Stan vb")
        return(rstan::vb(model, data = data, ...))
    }
}
