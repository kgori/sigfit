#' Fit mutational signatures
#'
#' \code{fit_signatures} performs MCMC sampling to fit a set of mutational signatures to a
#' collection of mutational catalogues and estimate the exposure of each catalogue to each signature.
#' Four models of signatures are available: multinomial, Poisson, normal and negative binomial. The
#' normal model can be used when \code{counts} contains continuous (non-integer) values, while the
#' negative binomial model is a more noise-robust version of the Poisson model.
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type.
#' @param signatures Mutational signatures to be fitted; either a numeric matrix with one row per
#' signature and one column per mutation type, or a list of matrices generated via
#' \code{\link{retrieve_pars}}.
#' @param exp_prior Numeric matrix with one row per sample and one column per signature, to be used
#' as the Dirichlet priors for the signature exposures. Default priors are uniform.
#' @param model Name of the model to sample from. Admits character values \code{"multinomial"}
#' (default), \code{"poisson"}, \code{"negbin"}, \code{"normal"}, \code{"nmf"} (an alias for
#' \code{"multinomial"}), and \code{"emu"} (an alias for \code{"poisson"}).
#' @param opportunities Numeric matrix of optional mutational opportunities, with one row per sample
#' and one column per mutation type. It also admits character values \code{"human-genome"} or
#' \code{"human-exome"}, in which case the mutational opportunities of the reference human
#' genome/exome will be used for every sample.
#' @param ... Additional arguments to be passed to \code{\link{rstan::sampling}}.
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{`data`}: list containing the input data supplied to the model.}
#'  \item{\code{`result`}: object of class stanfit, containing the output MCMC samples,
#'  as well as information about the model and the sampling process.}}
#' The model parameters (such as signature exposures) can be extracted from this
#' object using \code{\link{retrieve_pars}}.
#' @examples
#' \dontrun{
#' # Load example mutational catalogues and COSMIC signatures
#' data("counts_21breast")
#' data("cosmic_signatures_v2")
#'
#' # Fit signatures 1 to 4, using a custom prior that favors signature 1 over the rest
#' # (4 chains, 300 warmup iterations + 300 sampling iterations - use more in practice)
#' samples_1 <- fit_signatures(counts_21breast, cosmic_signatures_v2[1:4, ],
#'                             exp_prior = c(10, 1, 1, 1), iter = 600)
#'
#' # Fit all the signatures, running a single chain for many iterations
#' # (3000 warmup iterations + 10000 sampling iterations)
#' samples_2 <- fit_signatures(counts_21breast, cosmic_signatures_v2, chains = 1,
#'                             iter = 13000, warmup = 3000)
#' }
#' @importFrom "rstan" sampling
#' @export
fit_signatures <- function(counts, signatures, exp_prior = NULL, model = "multinomial",
                           opportunities = NULL, ...) {
    
    # Check argument restrictions
    if (!(model %in% c("normal", "poisson", "emu", "multinomial", "nmf", "negbin"))) {
        stop("'model' only admits values \"multinomial\", \"poisson\", \"normal\", \"negbin\", \"nmf\" and \"emu\".")
    }
    
    # Force counts and signatures to matrix
    counts_real <- to_matrix(counts)
    counts_int <- to_matrix(counts, int = TRUE)
    signatures <- to_matrix(signatures)
    
    if (model != "normal" & !all(counts_real - counts_int == 0)) {
        warning("Using a discrete model for non-integer counts: counts have been rounded.")
    }
    if (model == "normal" & all(counts_real - counts_int == 0)) {
        warning("Using a continuous model for integer counts; we recommend using a discrete model instead.")
    }
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    NSIG <- nrow(signatures)
    strand <- NCAT == 192  # strand bias indicator
    
    # Set up exposure priors and opportunities
    if (is.null(exp_prior)) {
        exp_prior <- matrix(1, nrow = NSAMP, ncol = NSIG)
    }
    opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)

    # Check dimensions are correct; should be
    # counts[NSAMP,NCAT], signatures[NSIG,NCAT], opportunities[NSAMP,NCAT], exp_prior[NSAMP,NSIG]
    stopifnot(ncol(signatures) == NCAT)
    stopifnot(all(dim(opportunities) == dim(counts)))
    stopifnot(all(dim(exp_prior) == c(NSAMP, NSIG)))

    # Select the model
    model_choices <- c("normal", "poisson", "emu", "multinomial", "nmf", "negbin")
    model <- match.arg(model, model_choices)
    if (model == "nmf") {
        cat("INFO: 'nmf' is an alias for 'multinomial'\n")
        model = "multinomial"
    }
    if (model == "emu") {
        cat("INFO: 'emu' is an alias for 'poisson'\n")
        model = "poisson"
    }

    # Construct model data list
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

#' Extract mutational signatures
#'
#' \code{extract_signatures} performs MCMC sampling to infer a set of mutational signatures
#' and their exposures from a collection of mutational catalogues. Four models of signatures are
#' available: multinomial, Poisson, normal and negative binomial. The normal model can be used when
#' \code{counts} contains continuous (non-integer) values, while the negative binomial model is a
#' more noise-robust version of the Poisson model. (However, the use of the negative binomial model
#' for signature extraction is discouraged due to its inefficiency.)
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type.
#' @param nsignatures Integer or integer vector indicating the number(s) of signatures to extract.
#' @param model Name of the model to sample from. Admits character values \code{"multinomial"}
#' (default), \code{"poisson"}, \code{"negbin"}, \code{"normal"}, \code{"nmf"} (an alias for
#' \code{"multinomial"}), and \code{"emu"} (an alias for \code{"poisson"}).
#' @param opportunities Numeric matrix of optional mutational opportunities, with one row per sample
#' and one column per mutation type. It also admits character values \code{"human-genome"} or
#' \code{"human-exome"}, in which case the mutational opportunities of the reference human
#' genome/exome will be used for every sample.
#' @param sig_prior Numeric matrix with one row per signature and one column per mutation type,
#' to be used as the Dirichlet priors for the mutational signatures. Only used when a single value
#' is provided for \code{nsignatures}. Default priors are uniform.
#' @param exp_prior Numeric matrix with one row per sample and one column per signature, to be used
#' as the Dirichlet priors for the signature exposures. Default priors are uniform.
#' @param dpp Logical indicating whether to use a Dirichlet process prior to infer the number of
#' mutational signatures (default is \code{FALSE}).
#' @param dpp_conc Numeric indicating the value of the concentration parameter for the Dirichlet
#' process prior (default is 1). Only used if \code{dpp=TRUE}.
#' @param stanfunc Character indicating the choice of rstan inference strategy.
#' Admits values \code{"sampling"}, \code{"optimizing"} and \code{"vb"}. The default value is
#' \code{"sampling"}, which corresponds to the full Bayesian MCMC approach. Alternatively,
#' \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates via numerical
#' optimization, while \code{"vb"} uses Variational Bayes to approximate the full posterior.
#' @param chains Integer indicating the number of chains used for MCMC (default is 1). The use of
#' multiple chains for signature extraction is discouraged, as it can result in an inference problem
#' called 'label switching'. This value is passed to \code{\link{rstan::sampling}}.
#' @param ... Additional arguments to be passed to the sampling function (by default,
#' \code{\link{rstan::sampling}}).
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{`data`}: list containing the input data supplied to the model.}
#'  \item{\code{`result`}: object of class stanfit, containing the output MCMC samples,
#'  as well as information about the model and the sampling process.}}
#' The model parameters (such as signatures and exposures) can be extracted from this object using
#' \code{\link{retrieve_pars}}.
#' If a range of numbers of signatures is provided via the \code{nsignatures} argument, a list is
#' returned in which the N-th element contains the extraction results for N signatures, as a list
#' with the structure described above.
#' @examples
#' \dontrun{
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
#' }
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @importFrom "rstan" extract
#' @export
extract_signatures <- function(counts, nsignatures, model = "multinomial", opportunities = NULL,
                               sig_prior = NULL, exp_prior = NULL, dpp = FALSE, dpp_conc = 1,
                               stanfunc = "sampling", chains = 1, ...) {
    
    # Check argument restrictions
    if (length(nsignatures) > 1 &
        (!is.null(sig_prior) | !is.null(exp_prior))) {
        stop("'sig_prior' and 'exp_prior' are only admitted when 'nsignatures' is a single integer.")
    }
    if (!(model %in% c("normal", "poisson", "emu", "multinomial", "nmf", "negbin"))) {
        stop("'model' only admits values \"multinomial\", \"poisson\", \"normal\", \"negbin\", \"nmf\" and \"emu\".")
    }
    if (length(chains) > 1 | !all.equal(chains, as.integer(chains))) {
        stop("'chains' must be a single positive integer.")
    }
    if (chains > 1) {
        warning("Using multiple chains for signature extraction can cause label switching problems and is discouraged.")
    }
    
    # Force counts to matrix
    counts_real <- to_matrix(counts)
    counts_int <- to_matrix(counts, int = TRUE)
    
    if (model != "normal" & !all(counts_real - counts_int == 0)) {
        warning("Using a discrete model for non-integer counts: counts have been rounded.")
    }
    if (model == "normal" & all(counts_real - counts_int == 0)) {
        warning("Using a continuous model for integer counts; we recommend using a discrete model instead.")
    }
    
    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    NSIG <- nsignatures[1]
    strand <- NCAT == 192  # strand bias indicator
    
    # Set up exposure priors and opportunities
    if (is.null(exp_prior)) {
        exp_prior <- matrix(1, nrow = NSAMP, ncol = NSIG)
    }
    if (is.null(sig_prior)) {
        sig_prior <- matrix(1, nrow = NSIG, ncol = NCAT)
    }
    opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)
    
    # Check dimensions are correct; should be
    # counts[NSAMP,NCAT], opportunities[NSAMP,NCAT], exp_prior[NSAMP,NSIG], sig_prior[NSIG,NCAT]
    stopifnot(all(dim(opportunities) == dim(counts)))
    stopifnot(all(dim(exp_prior) == c(NSAMP, NSIG)))
    stopifnot(all(dim(sig_prior) == c(NSIG, NCAT)))
    
    # Select the model
    model_choices <- c("normal", "poisson", "emu", "multinomial", "nmf", "negbin")
    model <- match.arg(model, model_choices)
    if (model == "nmf") {
        cat("INFO: 'nmf' is an alias for 'multinomial'\n")
        model = "multinomial"
    }
    if (model == "emu") {
        cat("INFO: 'emu' is an alias for 'poisson'\n")
        model = "poisson"
    }
    
    # Construct model data list
    dat <- list(
        C = NCAT,
        S = NSIG,
        G = NSAMP,
        counts = counts,
        counts_int = counts_int,
        counts_real = counts_real,
        opportunities = opportunities,
        kappa = exp_prior,
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
            dat$kappa <- matrix(1, nrow = NSAMP, ncol = n)
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
#' \code{fit_extract_signatures} performs MCMC sampling to simultaneously fit a set of 'fixed'
#' signatures to a collection of mutational catalogues (as in \code{\link{fit_signatures}}) and
#' extract a number of 'additional' signatures from the catalogues (as in
#' \code{\link{extract_signatures}}). Four models of signatures are available: multinomial, Poisson,
#' normal and negative binomial. The normal model can be used when \code{counts} contains continuous
#' (non-integer) values, while the negative binomial model is a more noise-robust version of the
#' Poisson model. (However, the use of the negative binomial model for signature extraction is
#' discouraged due to its inefficiency.)
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type.
#' @param signatures 'Fixed' mutational signatures to be fitted; either a numeric matrix with one
#' row per signature and one column per mutation type, or a list of matrices generated via
#' \code{\link{retrieve_pars}}.
#' @param num_extra_sigs Numeric indicating the number of 'additional' signatures to be extracted.
#' @param model Name of the model to sample from. Admits character values \code{"multinomial"}
#' (default), \code{"poisson"}, \code{"negbin"}, \code{"normal"}, \code{"nmf"} (an alias for
#' \code{"multinomial"}), and \code{"emu"} (an alias for \code{"poisson"}).
#' @param opportunities Numeric matrix of optional mutational opportunities, with one row per sample
#' and one column per mutation type. It also admits character values \code{"human-genome"} or
#' \code{"human-exome"}, in which case the mutational opportunities of the reference human
#' genome/exome will be used for every sample.
#' @param sig_prior Numeric matrix with one row per 'additional' signature and one column per
#' mutation type, to be used as the Dirichlet priors for the additional signatures to be extracted.
#' Default priors are uniform.
#' @param exp_prior Numeric matrix with one row per sample and one column per signature (including
#' both 'fixed' and 'additional' signatures), to be used as the Dirichlet priors for the signature
#' exposures. Default priors are uniform.
#' @param dpp Logical indicating whether to use a Dirichlet process prior to infer the number of
#' mutational signatures (default is \code{FALSE}).
#' @param dpp_conc Numeric indicating the value of the concentration parameter for the Dirichlet
#' process prior (default is 1). Only used if \code{dpp=TRUE}.
#' @param stanfunc Character indicating the choice of rstan inference strategy.
#' Admits values \code{"sampling"}, \code{"optimizing"} and \code{"vb"}. The default value is
#' \code{"sampling"}, which corresponds to the full Bayesian MCMC approach. Alternatively,
#' \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates via numerical
#' optimization, while \code{"vb"} uses Variational Bayes to approximate the full posterior.
#' @param chains Integer indicating the number of chains used for MCMC (default is 1). The use of
#' multiple chains for signature extraction is discouraged, as it can result in an inference problem
#' called 'label switching'. This value is passed to \code{\link{rstan::sampling}}.
#' @param ... Additional arguments to be passed to the sampling function (by default,
#' \code{\link{rstan::sampling}}).
#' @return A list with two elements:
#' \itemize{
#'  \item{\code{$data}: list containing the input data supplied to the model.}
#'  \item{\code{$result}: object of class stanfit, containing the output MCMC samples,
#'  as well as information about the model and sampling process.}}
#' The model parameters (such as signatures and exposures) can be extracted from this
#' object using \code{\link{retrieve_pars}}.
#' @examples
#' \dontrun{
#' # Simulate two catalogues using signatures 1, 4, 5, 7, with
#' # proportions 4:2:3:1 and 2:3:4:1, respectively
#' data("cosmic_signatures_v2")
#' probs <- rbind(c(0.4, 0.2, 0.3, 0.1) %*% cosmic_signatures_v2[c(1, 4, 5, 7), ],
#'                c(0.2, 0.3, 0.4, 0.1) %*% cosmic_signatures_v2[c(1, 4, 5, 7), ])
#' mutations <- rbind(t(rmultinom(1, 20000, probs[1, ])),
#'                    t(rmultinom(1, 20000, probs[2, ])))
#'
#' # Assuming that we do not know signature 7 a priori, but we know the others
#' # to be present, extract 1 signature while fitting signatures 1, 4 and 5.
#' # (400 warmup iterations + 400 sampling iterations - use more in practice)
#' mcmc_samples <- fit_extract_signatures(mutations, cosmic_signatures_v2[c(1, 4, 5), ],
#'                                        num_extra_sigs = 1, model = "nmf", iter = 800)
#'
#' # Plot original and extracted signature 7
#' extr_sigs <- retrieve_pars(mcmc_samples, "signatures")
#' plot_spectrum(cosmic_signatures_v2[7, ], pdf_path = "COSMIC_Sig7.pdf", name="COSMIC sig. 7")
#' plot_spectrum(extr_sigs, pdf_path = "Extracted_Sigs.pdf")
#' }
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @export
fit_extract_signatures <- function(counts, signatures, num_extra_sigs,
                                   model = "multinomial", opportunities = NULL, sig_prior = NULL,
                                   exp_prior = NULL, dpp = FALSE, dpp_conc = 1, 
                                   stanfunc = "sampling", chains = 1, ...) {

    # Check argument restrictions
    if (length(num_extra_sigs) != 1) {
        stop("'num_extra_sigs' must be a single positive integer.")
    }
    if (!(model %in% c("normal", "poisson", "emu", "multinomial", "nmf", "negbin"))) {
        stop("'model' only admits values \"multinomial\", \"poisson\", \"normal\", \"negbin\", \"nmf\" and \"emu\".")
    }
    if (length(chains) > 1 | !all.equal(chains, as.integer(chains))) {
        stop("'chains' must be a single positive integer.")
    }
    if (chains > 1) {
        warning("Using multiple chains for signature extraction can cause label switching problems and is discouraged.")
    }
    
    # Force counts and signatures to matrix
    counts_real <- to_matrix(counts)
    counts_int <- to_matrix(counts, int = TRUE)
    signatures <- to_matrix(signatures)
    
    if (model != "normal" & !all(counts_real - counts_int == 0)) {
        warning("Using a discrete model for non-integer counts: counts have been rounded.")
    }
    if (model == "normal" & all(counts_real - counts_int == 0)) {
        warning("Using a continuous model for integer counts; we recommend using a discrete model instead.")
    }
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    NSAMP <- nrow(counts)
    NCAT <- ncol(counts)
    NSIG <- nrow(signatures)
    strand <- NCAT == 192  # strand bias indicator
    
    # Set up exposure priors and opportunities
    if (is.null(exp_prior)) {
        exp_prior <- matrix(1, nrow = NSAMP, ncol = NSIG + num_extra_sigs)
    }
    if (is.null(sig_prior)) {
        sig_prior <- matrix(1, nrow = num_extra_sigs, ncol = NCAT)
    }
    opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)
    
    # Check dimensions are correct
    stopifnot(ncol(signatures) == NCAT)
    stopifnot(all(dim(opportunities) == dim(counts)))
    stopifnot(all(dim(exp_prior) == c(NSAMP, NSIG + num_extra_sigs)))
    stopifnot(all(dim(sig_prior) == c(num_extra_sigs, NCAT)))
    
    # Select the model
    model_choices <- c("normal", "poisson", "emu", "multinomial", "nmf", "negbin")
    model <- match.arg(model, model_choices)
    if (model == "nmf") {
        cat("INFO: 'nmf' is an alias for 'multinomial'\n")
        model = "multinomial"
    }
    if (model == "emu") {
        cat("INFO: 'emu' is an alias for 'poisson'\n")
        model = "poisson"
    }
    
    # Construct model data list
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
        kappa = exp_prior,
        concentration = dpp_conc,
        dpp = ifelse(dpp, 1, 0),
        family = switch(model,
                        nmf = 1, multinomial = 1, emu = 2, poisson = 2, negbin = 3, normal = 4)
    )
    
    cat("---\nFit-Ext: Fitting", NSIG, "signatures and extracting", num_extra_sigs,
        "signature(s) using", model, "model\n---\n")
    if (stanfunc == "sampling") {
        cat("Stan sampling:")
        out <- sampling(stanmodels$sigfit_fitext, data = dat, chains = chains,
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

#' Use optimization to generate initial parameter values for MCMC sampling
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type.
#' @param nsignatures Integer or integer vector indicating the number(s) of signatures to extract.
#' @param model Name of the model to sample from. Admits character values \code{"multinomial"}
#' (default), \code{"poisson"}, \code{"negbin"}, \code{"normal"}, \code{"nmf"} (an alias for
#' \code{"multinomial"}), and \code{"emu"} (an alias for \code{"poisson"}).
#' @param opportunities Numeric matrix of optional mutational opportunities, with one row per sample
#' and one column per mutation type. It also admits character values \code{"human-genome"} or
#' \code{"human-exome"}, in which case the mutational opportunities of the reference human
#' genome/exome will be used for every sample.
#' @param sig_prior Numeric matrix with one row per signature and one column per mutation type,
#' to be used as the Dirichlet priors for the mutational signatures. Only used when a single value
#' is provided for \code{nsignatures}. Default priors are uniform.
#' @param chains Integer indicating number of chains to be initialised (default is 1).
#' @param ... Additional arguments to be passed to \code{\link{rstan::optimizing}}.
#' @return List of initial values to be passed to \code{\link{extract_signatures}} via the
#' \code{init} argument.
#' @export
extract_signatures_initialiser <- function(counts, nsignatures, model = "multinomial",
                                           opportunities = NULL, sig_prior = NULL, exp_prior = NULL,
                                           chains = 1, ...) {

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
#' \code{get_initializer_list} extracts parameter values from an extraction which can be
#' used to initialise a follow-up extraction. One use case is to use a short, single-chain
#' extraction to initialise a longer, multi-chain extraction. This can mitigate against
#' label switching using the "Initialization around a single mode" strategy described in
#' the Stan documentation
#' \url{https://mc-stan.org/docs/2_18/stan-users-guide/label-switching-problematic-section.html}
#'
#' @param fitobj Result of a sigfit signature extraction.
#' @param chains Integer indicating the number of copies of the parameter list to be returned. Used
#' to initialize multiple chains.
#' @return A list of parameter lists.
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
