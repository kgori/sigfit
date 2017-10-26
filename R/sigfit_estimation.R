#' Fit mutational signatures
#' 
#' \code{fit_signatures} runs MCMC sampling to fit a set of mutational signatures 
#' to a collection of mutational catalogues and estimate the exposure of each
#' catalogue to each signature.
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param signatures Mutational signatures to be fitted. Either a matrix with one row per signature
#' and one column for each of the 96 mutation types, or a list of signatures generated via
#' \code{\link{retrieve_pars}}.
#' @param exp_prior Vector with one element per signature, to be used as the Dirichlet prior for 
#' the signature exposures in the sampling chain. Default prior is uniform (uninformative).
#' @param method Model to sample from; either \code{"nmf"} or \code{"emu"}.
#' @param opportunities Optional matrix of mutational opportunities for the "EMu" model 
#' (\code{method = "emu"}) method. Must be a matrix with same dimension as \code{counts}. 
#' If equal to \code{"human-genome"} or \code{"human-exome"}, the reference human genome/exome 
#' opportunities will be used for every sample.
#' @param cores Positive integer, number of cores to use. The default (\code{cores="auto"}) assigns
#' as many cores as detected in the system by \code{\link{parallel::detectCores()}}.
#' @param ... Arguments to pass to \code{rstan::sampling}.
#' @return A stanfit object containing the Monte Carlo samples from MCMC (from which the model
#' parameters can be extracted using \code{\link{retrieve_parameters}}), as well as information about
#' the model and sampling process.
#' @examples
#' # Load example mutational catalogues
#' data("counts_21breast")
#' 
#' # Fetch COSMIC signatures 
#' signatures <- fetch_cosmic_data()
#' 
#' # Fit signatures 1 to 4, using a custom prior that favors signature 1 over the rest
#' # (4 chains, 300 warmup iterations + 300 sampling iterations -- use more in practice)
#' samples_1 <- fit_signatures(counts_21breast, signatures[1:4, ], 
#'                             exp_prior = c(10, 1, 1, 1), iter = 600)
#' 
#' # Fit all the signatures, running a single chain for many iterations
#' # (3000 warmup iterations + 10000 sampling iterations)
#' samples_2 <- fit_signatures(counts_21breast, signatures, chains = 1, 
#'                             iter = 13000, warmup = 3000)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @export
fit_signatures <- function(counts, signatures, exp_prior = NULL, method = "nmf",
                           opportunities = NULL, cores = "auto", ...) {
    
    # Force counts and signatures to matrix
    counts <- to_matrix(counts)
    signatures <- to_matrix(signatures)
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    NSIG <- nrow(signatures)
    strand <- NCAT == 192  # strand bias indicator
    
    # Check exposure priors
    if (is.null(exp_prior)) {
        exp_prior = rep(1, NSIG)
    }
    exp_prior <- as.numeric(exp_prior)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NSIG, NCAT]
    stopifnot(ncol(signatures) == NCAT)
    stopifnot(length(exp_prior) == NSIG)
    
    if (method == "emu") {
        if (is.null(opportunities)) {
            opportunities <- matrix(1, nrow = NSAMP, ncol = NCAT)
        }
        else if (opportunities[1] == "human-genome") {
            opportunities <- build_opps_matrix(NSAMP, "genome", strand)
        }
        else if (opportunities[1] == "human-exome") {
            opportunities <- build_opps_matrix(NSAMP, "exome", strand)
        }
        if (!is.matrix(opportunities)) {
            opportunities <- as.matrix(opportunities)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        dat <- list(
            C = NCAT,
            S = NSIG,
            G = NSAMP,
            signatures = signatures,
            counts = counts,
            opps = opportunities,
            alpha = exp_prior
        )
        model <- stanmodels$sigfit_fit_emu
    }
    
    else {
        dat <- list(
            C = NCAT,
            S = NSIG,
            G = NSAMP,
            signatures = signatures,
            counts = counts,
            alpha = exp_prior
        )
        model <- stanmodels$sigfit_fit_nmf
    }
    
    # Set number of cores
    if (cores == "auto") {
        rstan_options(auto_write = TRUE)
        options(mc.cores = parallel::detectCores())
    }
    else if (cores %in% 2:parallel::detectCores()) {
        rstan_options(auto_write = TRUE)
        options(mc.cores = cores)
    }
    else if (!(cores %in% 1:parallel::detectCores())) {
        stop(paste0("'cores' must be a positive integer no larger than ", parallel::detectCores(), " (the number of cores in this system)."))
    }

    sampling(model, data = dat, pars = "multiplier", include = FALSE, ...)
}


#' Use optimization to generate initial parameter values for MCMC sampling
#' @export
extract_signatures_initialiser <- function(counts, nsignatures, method = "emu", opportunities = NULL, 
                                           sig_prior = NULL, chains = 1, ...) {
    
    opt <- extract_signatures(counts, nsignatures, method, opportunities, 
                              sig_prior, "optimizing", FALSE, ...)
    
    if (is.null(opt)) {
        warning("Parameter optimization failed - using random initialization")
        inits <- "random"
    }
    
    else {
        params <- list(
            signatures = matrix(opt$par[grepl("signatures", names(opt$par))], nrow = nsignatures),
            exposures_raw = matrix(opt$par[grepl("exposures_raw", names(opt$par))], nrow = nrow(counts))
        )
        inits = list()
        for (i in 1:chains) inits[[i]] <- p
    }
    inits
}

#' Extract signatures from a set of mutation counts
#' 
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param nsignatures Number (or vector of numbers) of signatures to extract.
#' @param method Either \code{"emu"} (default) or \code{"nmf"}.
#' @param opportunities Optional matrix of mutational opportunities for the "EMu" model 
#' (\code{method = "emu"}) method. Must be a matrix with same dimension as \code{counts}. 
#' If equal to \code{"human-genome"} or \code{"human-exome"}, the reference human genome/exome 
#' opportunities will be used for every sample.
#' @param sig_prior Matrix with one row per signature and one column per category, to be used as the Dirichlet 
#' priors for the signatures to be extracted. Only used when \code{nsignatures} is a scalar.
#' Default priors are uniform (uninformative).
#' @param stanfunc Choice of rstan inference strategy; admits values \code{"sampling"}, \code{"optimizing"}
#' and \code{"vb"}. \code{"sampling"} is the full Bayesian MCMC approach, and is the default. 
#' \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates via numerical optimization.
#' \code{"vb"} uses Variational Bayes to approximate the full posterior.
#' @param cores Positive integer, number of cores to use. The default (\code{cores="auto"}) assigns
#' as many cores as detected in the system by \code{\link{parallel::detectCores()}}.
#' @param ... Any other parameters to pass to the sampling function (by default, \code{\link{rstan::sampling}}).
#' (The number of chains is set to 1 and cannot be changed, to prevent 'label switching' problems.)
#' @return A stanfit object containing the Monte Carlo samples from MCMC (from which the model
#' parameters can be extracted using \code{\link{retrieve_parameters}}), as well as information about
#' the model and sampling process.
#' @examples
#' # Load example mutational catalogues
#' data("counts_21breast")
#' 
#' # Extract 2 to 6 signatures using the NMF (multinomial) model
#' # (400 warmup iterations + 400 sampling iterations -- use more in practice)
#' samples_nmf <- extract_signatures(counts_21breast, nsignatures = 2:6, 
#'                                   method = "nmf", iter = 800)
#' 
#' # Extract 4 signatures using the EMu (Poisson) model
#' # (400 warmup iterations + 800 sampling iterations -- use more in practice)
#' samples_emu <- extract_signatures(counts_21breast, nsignatures = 4, method = "emu", 
#'                                   opportunities = "human-genome",
#'                                   iter = 1200, warmup = 400)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @importFrom "rstan" extract
#' @export
extract_signatures <- function(counts, nsignatures, method = "emu", opportunities = NULL, 
                               sig_prior = NULL, stanfunc = "sampling", cores = "auto", ...) {
    
    if (!is.null(sig_prior) & length(nsignatures) > 1) {
        stop("'sig_prior' is only admitted when 'nsignatures' is a scalar (single value).")
    }
    
    # Force counts to matrix
    counts <- to_matrix(counts)
    
    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    strand <- NCAT == 192  # strand bias indicator
    
    # EMu model
    if (method == "emu") {
        # Build opportunities matrix
        if (is.null(opportunities)) {
            opportunities <- matrix(1, nrow = NSAMP, ncol = NCAT)
        }
        else if (opportunities[1] == "human-genome") {
            opportunities <- build_opps_matrix(NSAMP, "genome", strand)
        }
        else if (opportunities[1] == "human-exome") {
            opportunities <- build_opps_matrix(NSAMP, "exome", strand)
        }
        if (!is.matrix(opportunities)) {
            opportunities <- as.matrix(opportunities)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        model <- stanmodels$sigfit_ext_emu
        
        dat <- list(
            C = NCAT,
            G = NSAMP,
            S = as.integer(nsignatures[1]),
            counts = counts,
            opps = opportunities,
            alpha = sig_prior
        )
    }
    
    # NMF model
    else if (method == "nmf") {
        if (!is.null(opportunities)) {
            warning("Using \"nmf\" model: 'opportunities' will be ignored.")
        }
        
        model <- stanmodels$sigfit_ext_nmf
        
        dat <- list(
            C = NCAT,
            G = NSAMP,
            S = as.integer(nsignatures[1]),
            counts = counts,
            alpha = sig_prior
        )
    }
    else {
        stop("'method' must be either \"emu\" or \"nmf\".")
    }
    
    # Set number of cores
    if (cores == "auto") {
        rstan_options(auto_write = TRUE)
        options(mc.cores = parallel::detectCores())
    }
    else if (cores %in% 2:parallel::detectCores()) {
        rstan_options(auto_write = TRUE)
        options(mc.cores = cores)
    }
    else if (!(cores %in% 1:parallel::detectCores())) {
        stop(paste0("'cores' must be a positive integer no larger than ", parallel::detectCores(), " (the number of cores in this system)."))
    }
    
    # Extract signatures for each nsignatures value
    if (length(nsignatures) > 1) {
        out <- vector(mode = "list", length = max(nsignatures))
        for (n in nsignatures) {
            # Complete model data
            dat$alpha <- matrix(1, nrow = n, ncol = NCAT)
            dat$S <- as.integer(n)
            
            cat("Extracting", n, "signatures\n")
            if (stanfunc == "sampling") {
                cat("Stan sampling:")
                out[[n]] <- sampling(model, data = dat, chains = 1, ...)
                
            }
            else if (stanfunc == "optimizing") {
                cat("Stan optimizing:")
                out[[n]] <- optimizing(model, data = dat, ...)
            }
            else if (stanfunc == "vb") {
                cat("Stan vb:")
                out[[n]] <- vb(model, data = dat, ...)
            }
        }
        
        names(out) <- paste0("nsignatures=", 1:length(out))
        
        # Plot goodness of fit and best number of signatures
        out$best <- plot_gof(out, counts)
    }
    
    # Single nsignatures value case
    else {
        # Check signature priors
        if (is.null(sig_prior)) {
            sig_prior <- matrix(1, nrow = nsignatures, ncol = NCAT)
        }
        sig_prior <- as.matrix(sig_prior)
        stopifnot(nrow(sig_prior) == nsignatures)
        stopifnot(ncol(sig_prior) == NCAT)
        dat$alpha <- sig_prior
        
        cat("Extracting", nsignatures, "signatures\n")
        if (stanfunc == "sampling") {
            cat("Stan sampling:")
            out <- sampling(model, data = dat, chains = 1, ...)
        }
        else if (stanfunc == "optimizing") {
            cat("Stan optimizing:")
            out <- optimizing(model, data = dat, ...)
        }
        else if (stanfunc == "vb") {
            cat("Stan vb:")
            out <- vb(model, data = dat, ...)
        }
    }
    out
}

#' Fit-and-extract mutational signatures
#' 
#' \code{fit_extract_signatures} fits signatures to estimate exposures in a set of mutation counts
#' and extracts additional signatures present in the samples.
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param signatures Fixed mutational signatures (columns) to be fitted. Either a matrix with one row 
#' per signature and one column for each of the 96 mutation types, or a list of signatures generated via
#' \code{\link{retrieve_pars}}.
#' @param num_extra_sigs Number of additional signatures to be extracted.
#' @param sig_prior Matrix with one row per additional signature and one column per category, to be used as the
#' Dirichlet priors for the additional signatures to be extracted. Default priors are uniform (uninformative).
#' @param stanfunc \code{"sampling"}|\code{"optimizing"}|\code{"vb"} Choice of rstan 
#' inference strategy. \code{"sampling"} is the full Bayesian MCMC approach, and is the 
#' default. \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates 
#' via numerical optimization. \code{"vb"} uses Variational Bayes to approximate the 
#' full posterior.
#' @param cores Positive integer, number of cores to use. The default (\code{cores="auto"}) assigns
#' as many cores as detected in the system by \code{\link{parallel::detectCores()}}.
#' @param ... Any other parameters to pass to the sampling function (by default, \code{\link{rstan::sampling}}).
#' (The number of chains is set to 1 and cannot be changed, to prevent 'label switching' problems.)
#' @return A stanfit object containing the Monte Carlo samples from MCMC (from which the model
#' parameters can be extracted using \code{\link{retrieve_parameters}}), as well as information about
#' the model and sampling process.
#' @examples
#' # Fetch COSMIC signatures 
#' signatures <- fetch_cosmic_data()
#'
#' # Simulate two catalogues using signatures 1, 4, 5, 7, with
#' # proportions 4:2:3:1 and 2:3:4:1, respectively
#' probs <- rbind(c(0.4, 0.2, 0.3, 0.1) %*% signatures[c(1, 4, 5, 7), ],
#'                c(0.2, 0.3, 0.4, 0.1) %*% signatures[c(1, 4, 5, 7), ])
#' mutations <- rbind(t(rmultinom(1, 20000, probs[1, ])),
#'                    t(rmultinom(1, 20000, probs[2, ])))  
#' 
#' # Assuming that we do not know signature 7 a priori, but we know the others
#' # to be present, extract 1 signature while fitting signatures 1, 4 and 5.
#' # (400 warmup iterations + 400 sampling iterations -- use more in practice)
#' mcmc_samples <- fit_extract_signatures(mutations, signatures = signatures[c(1, 4, 5), ],
#'                                        num_extra_sigs = 1, method = "nmf", iter = 800)
#' 
#' # Plot original and extracted signature 7
#' extr_sigs <- retrieve_pars(mcmc_samples, "signatures")
#' plot_spectrum(signatures[7, ], pdf_path = "COSMIC_Sig7.pdf", name="COSMIC sig. 7")
#' plot_spectrum(extr_sigs, pdf_path = "Extracted_Sigs.pdf")
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @export
fit_extract_signatures <- function(counts, signatures, num_extra_sigs, 
                                   method = "nmf", opportunities = NULL, sig_prior = NULL,
                                   stanfunc = "sampling", cores = "auto", ...) {
    # Check that num_extra_sigs is scalar
    if (length(num_extra_sigs) > 1) {
        stop("'num_extra_sigs' must be an integer scalar.")
    }
    
    
    # Force counts and signatures to matrix
    counts <- to_matrix(counts)
    signatures <- to_matrix(signatures)
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    NSIG <- nrow(signatures)
    strand <- NCAT == 192  # strand bias indicator
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NSIG, NCAT]
    stopifnot(ncol(signatures) == NCAT)
    
    # Check signature priors
    if (is.null(sig_prior)) {
        sig_prior <- matrix(1, nrow = num_extra_sigs, ncol = NCAT)
    }
    sig_prior <- as.matrix(sig_prior)
    stopifnot(nrow(sig_prior) == num_extra_sigs)
    stopifnot(ncol(sig_prior) == NCAT)
    
    # EMu model
    if (method == "emu") {
        # Build opportunities matrix
        if (is.null(opportunities)) {
            opportunities <- matrix(1, nrow = NSAMP, ncol = NCAT)
        }
        else if (is.character(opportunities) & opportunities == "human-genome") {
            opportunities <- build_opps_matrix(NSAMP, "genome", strand)
        }
        else if (is.character(opportunities) & opportunities == "human-exome") {
            opportunities <- build_opps_matrix(NSAMP, "exome", strand)
        }
        if (!is.matrix(opportunities)) {
            opportunities <- as.matrix(opportunities)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        model <- stanmodels$sigfit_fitex_emu
        dat <- list(
            C = NCAT,
            S = NSIG,
            G = NSAMP,
            N = as.integer(num_extra_sigs),
            fixed_sigs = signatures,
            counts = counts,
            opps = opportunities,
            alpha = sig_prior
        )
    }
    
    # NMF model
    else if (method == "nmf") {
        model <- stanmodels$sigfit_fitex_nmf
        dat <- list(
            C = NCAT,
            S = NSIG,
            G = NSAMP,
            N = as.integer(num_extra_sigs),
            fixed_sigs = signatures,
            counts = counts,
            alpha = sig_prior
        )
    }
    
    # Set number of cores
    if (cores == "auto") {
        rstan_options(auto_write = TRUE)
        options(mc.cores = parallel::detectCores())
    }
    else if (cores %in% 2:parallel::detectCores()) {
        rstan_options(auto_write = TRUE)
        options(mc.cores = cores)
    }
    else if (!(cores %in% 1:parallel::detectCores())) {
        stop(paste0("'cores' must be a positive integer no larger than ", parallel::detectCores(), " (the number of cores in this system)."))
    }
    
    if (stanfunc == "sampling") {
        cat("Stan sampling:")
        sampling(model, data = dat, chains = 1, 
                 pars = c("extra_sigs", "probs", "exposures_raw", "lambda"),
                 include = FALSE, ...)
    }
    else if (stanfunc == "optimizing") {
        cat("Stan optimizing:")
        optimizing(model, data = dat, ...)
    }
    else if (stanfunc == "vb") {
        cat("Stan vb")
        vb(model, data = dat, ...)
    }
}
