#' Fit mutational signatures
#'
#' \code{fit_signatures} runs MCMC sampling to fit a set of mutational signatures
#' to a collection of mutational catalogues and estimate the exposure of each
#' catalogue to each signature.
#' @param counts Integer matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Any decimal values will be rounded to integers.
#' @param signatures Mutational signatures to be fitted. Either a numeric matrix with one row per signature
#' and one column per mutation type, or a list of matrices generated via
#' \code{\link{retrieve_pars}}.
#' @param exp_prior Numeric vector with one element per signature, to be used as the Dirichlet prior for
#' the signature exposures in the sampling chain. Default prior is uniform (uninformative).
#' @param model Character; model to sample from. Admits values \code{"nmf"} (default) or \code{"emu"}.
#' @param opportunities Numeric matrix of optional mutational opportunities for the "EMu" model
#' (\code{model = "emu"}). Must be a matrix with same dimension as \code{counts}.
#' Alternatively, it also admits character values \code{"human-genome"} or \code{"human-exome"},
#' in which case the reference human genome/exome opportunities will be used for every sample.
#' @param ... Arguments to pass to \code{rstan::sampling}.
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{$data}: list containing the input data supplied to the model.}
#'  \item{\code{$result}: object of class stanfit, containing the output MCMC samples,
#'  as well as information about the model and sampling process.}}
#' The model parameters (such as signatures and exposures) can be extracted from this
#' object using \code{\link{retrieve_pars}}.
#' @examples
#' \dontrun{
#' # Load example mutational catalogues and COSMIC signatures
#' data("counts_21breast")
#' data("cosmic_signatures")
#'
#' # Fit signatures 1 to 4, using a custom prior that favors signature 1 over the rest
#' # (4 chains, 300 warmup iterations + 300 sampling iterations - use more in practice)
#' samples_1 <- fit_signatures(counts_21breast, cosmic_signatures[1:4, ],
#'                             exp_prior = c(10, 1, 1, 1), iter = 600)
#'
#' # Fit all the signatures, running a single chain for many iterations
#' # (3000 warmup iterations + 10000 sampling iterations)
#' samples_2 <- fit_signatures(counts_21breast, cosmic_signatures, chains = 1,
#'                             iter = 13000, warmup = 3000)
#' }
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @export
fit_signatures <- function(counts, signatures, exp_prior = NULL, model = "multinomial",
                           opportunities = NULL, ...) {
    
    # Force counts and signatures to matrix
    counts_real <- to_matrix(counts)
    counts_int <- to_matrix(counts, int = TRUE)
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

    # Select the model
    model_choices <- c("normal", "poisson", "emu", "multinomial", "nmf", "negbin")
    model <- match.arg(model, model_choices)
    if (model == "nmf") {
        cat("INFO:'nmf' is an alias for 'multinomial'\n")
        model = "multinomial"
    }
    if (model == "emu") {
        cat("INFO:'emu' is an alias for 'poisson'\n")
        model = "poisson"
    }

    # Set up the opportunities
    opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)
    stopifnot(all(dim(opportunities) == dim(counts)))

    dat <- list(
        C = NCAT,
        S = NSIG,
        G = NSAMP,
        signatures = signatures,
        counts_int = counts_int,
        counts_real = counts_real,
        kappa = exp_prior,
        opportunities = opportunities,
        family = switch(model,
                        nmf = 1, multinomial = 1, emu = 2, poisson = 2, negbin = 3, normal = 4)
    )

    cat("---\nFitting", NSIG, "signatures using", model, "model\n---\n")
    cat("Stan sampling:")
    out <- sampling(stanmodels$sigfit_fit, data = dat, pars = "multiplier", include = FALSE, ...)

    list("data" = dat,
         "result" = out)
}


#' Use optimization to generate initial parameter values for MCMC sampling
#' @param counts Integer matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Any decimal values will be rounded to integers.
#' @param nsignatures Integer or integer vector indicating the number(s) of signatures to extract.
#' @param model Character indicating the model to sample from. Admits values \code{"nmf"} or
#' \code{"emu"} (the default).
#' @param opportunities Numeric matrix of optional mutational opportunities for the "EMu" model
#' (\code{model = "emu"}). Must be a matrix with same dimension as \code{counts}.
#' Alternatively, it also admits character values \code{"human-genome"} or \code{"human-exome"},
#' in which case the reference human genome/exome opportunities will be used for every sample.
#' @param sig_prior Numeric matrix with one row per signature and one column per mutation type,
#' to be used as the Dirichlet priors for the signatures to be extracted. Only used when
#' \code{nsignatures} is a scalar. Default priors are uniform (uninformative).
#' @param chains Integer indicating number of chains to be initialised (default is 1).
#' @param ... Additional arguments to pass to \code{rstan::optimizing}.
#' @return List of initial values to be passed to \code{\link{extract_signatures}} via the
#' \code{init} argument.
#' @export
extract_signatures_initialiser <- function(counts, nsignatures, model = "emu", opportunities = NULL,
                                           sig_prior = NULL, exp_prior = 1, chains = 1, ...) {

    opt <- extract_signatures(counts, nsignatures, model, opportunities,
                              sig_prior, exp_prior, stanfunc = "optimizing", ...)

    if (is.null(opt)) {
        stop("Parameter optimization failed")
    }
    else {
        return (get_initializer_list(opt, chains))
    }
}

#' Query a signature extraction result for parameter values that can initialise a follow-up
#' extraction.
#'
#' \code{get_initializer_list} extracts parameter values from an extraction that can be
#' used to initialise a follow-up extraction. One use case is to use a short, single-chain
#' extraction to initialise a longer, multi-chain extraction. This can mitigate against
#' label switching using the "Initialization around a single mode" strategy described in
#' the Stan documentation
#' \url{https://mc-stan.org/docs/2_18/stan-users-guide/label-switching-problematic-section.html}
#'
#' @param fitobj Result of a sigfit signature extraction
#' @param chains (integer) return a copy of the parameter list for each of \code{chains}. Used
#' to initialize multiple chains.
#' @return list of parameter lists.
#' @export
get_initializer_list <- function(fitobj, chains = 1) {
    params <- list(
        activities = as.matrix(retrieve_pars(fitobj, "activities")$mean),
        exposures = as.matrix(retrieve_pars(fitobj, "exposures")$mean),
        signatures = as.matrix(retrieve_pars(fitobj, "signatures")$mean)
    )

    inits <- list()
    for (i in 1:chains) inits[[i]] <- params
    inits
}

#' Extract signatures from a set of mutation counts
#'
#' \code{extract_signatures} runs MCMC sampling to extract a set of mutational signatures
#' and their exposures from a collection of mutational catalogues.
#' @param counts Integer matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Any decimal values will be rounded to integers.
#' @param nsignatures Integer or integer vector indicating the number(s) of signatures to extract.
#' @param model Character indicating the model to sample from. Admits values \code{"nmf"}
#' (the default) or \code{"emu"}.
#' @param opportunities Numeric matrix of optional mutational opportunities for the "EMu" model
#' (\code{model = "emu"}). Must be a matrix with same dimension as \code{counts}.
#' Alternatively, it also admits character values \code{"human-genome"} or \code{"human-exome"},
#' in which case the reference human genome/exome opportunities will be used for every sample.
#' @param sig_prior Numeric matrix with one row per signature and one column per mutation type,
#' to be used as the Dirichlet priors for the signatures to be extracted. Only used when
#' \code{nsignatures} is a scalar. Default priors are uniform (uninformative).
#' @param exp_prior Numeric; hyperparameter of the Dirichlet prior given to the exposures.
#' Default value is 1 (uniform, uninformative).
#' @param stanfunc Character indicating the choice of rstan inference strategy.
#' Admits values \code{"sampling"}, \code{"optimizing"} and \code{"vb"}. The default value is
#' \code{"sampling"}, which corresponds to the full Bayesian MCMC approach. Alternatively,
#' \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates via numerical
#' optimization, while \code{"vb"} uses Variational Bayes to approximate the full posterior.
#' @param ... Any additional parameters to be passed to the sampling function (by default,
#' \code{rstan::sampling}). (Note that the number of chains is set to 1 and cannot
#' be changed, in order to prevent 'label switching' problems.)
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{$data}: list containing the input data supplied to the model.}
#'  \item{\code{$result}: object of class stanfit, containing the output MCMC samples,
#'  as well as information about the model and sampling process.}}
#' The model parameters (such as signatures and exposures) can be extracted from this
#' object using \code{\link{retrieve_pars}}.
#' If a range of numbers of signatures is provided via the \code{nsignatures} argument, a list is returned where
#' the N-th element contains the extraction results for N signatures, as a list with the structure described above.
#' @examples
#' # Load example mutational catalogues
#' data("counts_21breast")
#'
#' # Extract 2 to 6 signatures using the NMF (multinomial) model
#' # (400 warmup iterations + 400 sampling iterations - use more in practice)
#' samples_nmf <- extract_signatures(counts_21breast, nsignatures = 2:6,
#'                                   model = "nmf", iter = 800)
#'
#' # Extract 4 signatures using the EMu (Poisson) model
#' # (400 warmup iterations + 800 sampling iterations - use more in practice)
#' samples_emu <- extract_signatures(counts_21breast, nsignatures = 4, model = "emu",
#'                                   opportunities = "human-genome",
#'                                   iter = 1200, warmup = 400)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @importFrom "rstan" extract
#' @export
extract_signatures <- function(counts, nsignatures, model = "multinomial", opportunities = NULL,
                               sig_prior = NULL, exp_prior = 1, dpp = FALSE, dpp_conc = 1, stanfunc = "sampling",
                               chains = 1, ...) {
    if (!is.null(sig_prior) & length(nsignatures) > 1) {
        stop("'sig_prior' is only admitted when 'nsignatures' is a scalar (single value).")
    }

    if (length(chains) > 1 | !all.equal(chains, as.integer(chains))) {
        stop("'chains' must be a single integer value")
    }

    stopifnot(is.numeric(exp_prior) & exp_prior > 0)
    if (length(exp_prior) > 1 | exp_prior != 1) warning("Setting custom priors on exposures is not implemented yet")

    # Force counts to matrix
    counts_real <- to_matrix(counts)
    counts_int <- to_matrix(counts, int = TRUE)

    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    strand <- NCAT == 192  # strand bias indicator

    # Select the model
    model_choices <- c("normal", "poisson", "emu", "multinomial", "nmf", "negbin")
    model <- match.arg(model, model_choices)
    if (model == "nmf") {
        cat("INFO:'nmf' is an alias for 'multinomial'\n")
        model = "multinomial"
    }
    if (model == "emu") {
        cat("INFO:'emu' is an alias for 'poisson'\n")
        model = "poisson"
    }

    # Set up the opportunities
    opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)
    stopifnot(all(dim(opportunities) == dim(counts)))
    
    dat <- list(
        C = NCAT,
        S = nsignatures,
        G = NSAMP,
        counts = counts,
        counts_int = counts_int,
        counts_real = counts_real,
        kappa = exp_prior,
        opportunities = opportunities,
        alpha = sig_prior,
        family = switch(model,
                        nmf = 1, multinomial = 1, emu = 2, poisson = 2, negbin = 3, normal = 4),
        concentration = dpp_conc,
        dpp = ifelse(dpp, 1, 0)
    )

    # Extract signatures for each nsignatures value
    if (length(nsignatures) > 1) {
        out <- vector(mode = "list", length = max(nsignatures))
        for (n in nsignatures) {
            # Complete model data
            dat$alpha <- matrix(1, nrow = n, ncol = NCAT)
            dat$kappa <- rep(exp_prior, n)
            dat$S <- as.integer(n)

            cat("---\nExtracting", n, "signatures using", model, "model\n---\n")
            if (stanfunc == "sampling") {
                cat("Stan sampling:")
                out[[n]] <- list("data" = dat,
                                 "result" = sampling(stanmodels$sigfit_ext,
                                                     data = dat, chains = chains, ...))
            }
            else if (stanfunc == "optimizing") {
                cat("Stan optimizing:")
                out[[n]] <- list("data" = dat,
                                 "result" = optimizing(stanmodels$sigfit_ext,
                                                       data = dat, ...))
            }
            else if (stanfunc == "vb") {
                cat("Stan vb:")
                out[[n]] <- list("data" = dat,
                                 "result" = vb(stanmodels$sigfit_ext,
                                               data = dat, ...))
            }
        }

        names(out) <- paste0("nsignatures=", 1:length(out))

        # Plot goodness of fit and best number of signatures
        out$best <- plot_gof(out)
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
        dat$kappa <- rep(exp_prior, nsignatures)

        cat("---\nExtracting", nsignatures, "signatures using", model, "model\n---\n")
        if (stanfunc == "sampling") {
            cat("Stan sampling:")
            out <- list("data" = dat,
                        "result" = sampling(stanmodels$sigfit_ext, data = dat, chains = chains, ...))
        }
        else if (stanfunc == "optimizing") {
            cat("Stan optimizing:")
            out <- list("data" = dat,
                        "result" = optimizing(stanmodels$sigfit_ext,
                                              data = dat, ...))
        }
        else if (stanfunc == "vb") {
            cat("Stan vb:")
            out <- list("data" = dat,
                        "result" = vb(stanmodels$sigfit_ext,
                                      data = dat, ...))
        }
    }

    out
}

#' Fit-and-extract mutational signatures
#'
#' \code{fit_extract_signatures} fits signatures to estimate exposures in a set of mutation counts
#' and simultaneously extracts additional signatures present in the samples.
#' @param counts Integer matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Any decimal values will be rounded to integers.
#' @param signatures Fixed mutational signatures to be fitted. Either a numeric matrix with one row
#' per signature and one column per mutation type, or a list of matrices generated via
#' \code{\link{retrieve_pars}}.
#' @param num_extra_sigs Numeric indicating the number of additional signatures to be extracted.
#' @param model Character indicating the model to sample from. Admits values \code{"nmf"} (default) or \code{"emu"}.
#' @param opportunities Numeric matrix of optional mutational opportunities for the "EMu" model
#' (\code{model = "emu"}). Must be a matrix with same dimension as \code{counts}.
#' Alternatively, it also admits character values \code{"human-genome"} or \code{"human-exome"},
#' in which case the reference human genome/exome opportunities will be used for every sample.
#' @param sig_prior Numeric matrix with one row per additional signature and one column per category, to be used as the
#' Dirichlet priors for the additional signatures to be extracted. Default priors are uniform (uninformative).
#' @param exp_prior Numeric; hyperparameter of the Dirichlet prior given to the exposures.
#' Default value is 1 (uniform, uninformative).
#' @param stanfunc Character indicating the choice of rstan inference strategy. Admits values \code{"sampling"},
#' \code{"optimizing"} and \code{"vb"}. The default value is
#' \code{"sampling"}, which corresponds to the full Bayesian MCMC approach. Alternatively,
#' \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates via numerical
#' optimization, while \code{"vb"} uses Variational Bayes to approximate the full posterior.
#' @param ... Any additional parameters to pass to the sampling function (by default, \code{rstan::sampling}).
#' (Note that the number of chains is set to 1 and cannot be changed, to prevent 'label switching' problems.)
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{$data}: list containing the input data supplied to the model.}
#'  \item{\code{$result}: object of class stanfit, containing the output MCMC samples,
#'  as well as information about the model and sampling process.}}
#' The model parameters (such as signatures and exposures) can be extracted from this
#' object using \code{\link{retrieve_pars}}.
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
#' # (400 warmup iterations + 400 sampling iterations - use more in practice)
#' mcmc_samples <- fit_extract_signatures(mutations, signatures = signatures[c(1, 4, 5), ],
#'                                        num_extra_sigs = 1, model = "nmf", iter = 800)
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
                                   model = "multinomial", opportunities = NULL, sig_prior = NULL,
                                   exp_prior = 1, dpp = FALSE, dpp_conc = 1, stanfunc = "sampling", ...) {
    # Check that num_extra_sigs is scalar
    if (length(num_extra_sigs) != 1) {
        stop("'num_extra_sigs' must be an integer scalar.")
    }

    stopifnot(is.numeric(exp_prior) & exp_prior > 0)

    # Force counts and signatures to matrix
    counts_real <- to_matrix(counts)
    counts_int <- to_matrix(counts, int = TRUE)
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

    # Select the model
    model_choices <- c("normal", "poisson", "emu", "multinomial", "nmf", "negbin")
    model <- match.arg(model, model_choices)
    if (model == "nmf") {
        cat("INFO:'nmf' is an alias for 'multinomial'\n")
        model = "multinomial"
    }
    if (model == "emu") {
        cat("INFO:'emu' is an alias for 'poisson'\n")
        model = "poisson"
    }

    # Set up the opportunities
    opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)
    stopifnot(all(dim(opportunities) == dim(counts)))

    dat <- list(
        C = NCAT,
        S = NSIG,
        G = NSAMP,
        N = as.integer(num_extra_sigs),
        fixed_sigs = signatures,
        counts_int = counts_int,
        counts_real = counts_real,
        opportunities = opportunities,
        alpha = sig_prior,
        kappa = rep(exp_prior, NSIG + num_extra_sigs),
        concentration = dpp_conc,
        dpp = ifelse(dpp, 1, 0),
        family = switch(model,
                        nmf = 1, multinomial = 1, emu = 2, poisson = 2, negbin = 3, normal = 4)
    )
    
    cat("---\nFit-Ext: Fitting", NSIG, "signatures and extracting", num_extra_sigs,
        "signature(s) using", model, "model\n---\n")
    if (stanfunc == "sampling") {
        cat("Stan sampling:")
        out <- sampling(stanmodels$sigfit_fitext, data = dat, chains = 1,
                        pars = "extra_sigs", include = FALSE, ...)
    }
    else if (stanfunc == "optimizing") {
        cat("Stan optimizing:")
        out <- optimizing(stanmodels$sigfit_fitext, data = dat, ...)
    }
    else if (stanfunc == "vb") {
        cat("Stan vb")
        out <- vb(stanmodels$sigfit_fitext, data = dat, ...)
    }

    list("data" = dat,
         "result" = out)
}
