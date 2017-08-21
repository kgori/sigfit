#' Deal with zero values in a signature
remove_zeros_ <- function(mtx, min_allowed = 1e-9) {
    t(apply(mtx, 1, function(row) {
        row[row < min_allowed] <- min_allowed
        row / sum(row)
    }))
}

#' Reverse complement a nucleotide sequence string
#' Input is string, output is character vector
rev_comp <- function(nucleotides) {
    as.character(rev(sapply(strsplit(nucleotides, "")[[1]], function(nuc) {
        if (nuc == "A") "T"
        else if (nuc == "C") "G"
        else if (nuc == "G") "C"
        else if (nuc == "T") "A"
    })))
}

#' Returns character vector of 96 (pyrimidine) trinucleotide mutation types
mut_types <- function() {
    bases <- c("A", "C", "G", "T")
    paste0(rep(rep(bases, each = 4), 6),
           rep(bases[c(2, 4)], each = 48),
           rep(bases, 6 * 16 / 4),
           ">",
           rep(rep(bases, each = 4), 6),
           c(rep(bases[-2], each = 16), rep(bases[-4], each = 16)),
           rep(bases, 6 * 16 / 4))
}

#' Builds a mutational catalogue from a table containing the base change and
#' trinucleotide context of each single-nucleotide variant.
#' @param variants Table with one row per single-nucleotide variant, and three columns:
#' REF base (character in {A,C,G,T}); ALT base (character in {A,C,G,T}); and trinucleotide
#' context at the variant location (sequence between the positions immediately before and 
#' after the variant, in string format; e.g. "TCA").
#' @export
build_catalogue <- function(variants) {
    # Check that REF base coincides with middle base in trinucleotide
    if (any(variants[,1] != sapply(strsplit(variants[,3], split=""), function(x) x[2])))
        stop("REF base (first column) must always be equal to middle base of the trinucleotide context (third column).")
    
    # Obtain mutation types, collapsed such that they refer to pyrimidine bases
    vars_collapsed <- apply(variants, 1, function(var) {
        if (var[1] %in% c("C", "T")) {
            trinuc <- strsplit(var[3], "")[[1]]
            alt <- var[2]
        }
        else {
            trinuc <- rev_comp(var[3])
            alt <- rev_comp(var[2])
        }
        paste0(paste(trinuc, collapse=""), ">", trinuc[1], alt, trinuc[3])
    })

    # Count number of occurrences of each mutation type
    sapply(mut_types(), function(type) {
        sum(grepl(type, vars_collapsed, fixed = TRUE))
    })
}

#' Fetches COSMIC's estimated mutational signatures
#' @param reorder Reorders the matrix by substitution type and trinucleotide
#' @export
fetch_cosmic_data <- function(reorder = TRUE, remove_zeros = TRUE) {
    cosmic_sigs <- read.table('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt', 
                              header = TRUE, sep = '\t', check.names = FALSE)
    if (reorder) {
        cosmic_sigs <- cosmic_sigs[order(cosmic_sigs[["Substitution Type"]], cosmic_sigs[["Trinucleotide"]]),]
    }
    rownames(cosmic_sigs) <- cosmic_sigs[["Somatic Mutation Type"]]
    cosmic_sigs <- t(cosmic_sigs[, paste("Signature", 1:30)])
    if (remove_zeros) cosmic_sigs <- remove_zeros_(cosmic_sigs)
    cosmic_sigs
}

#' Returns human genome or exome trinucleotide frequencies. This is useful to
#' de-normalize signatures that were obtained using the mutational opportunities
#' from the human genome/exome, in order to compare them to signatures extracted
#' without incorporating opportunities.
#' @param type Either "genome" (default) or "exome".
#' @examples
#' Extract signatures using human exome opportunitites
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "emu", 
#' opportunities = "human-exome")
#' sigs <- retrieve_pars(samples, "signatures")
#' 
#' # De-normalize (mean) extracted signatures
#' freqs <- human_trinuc_freqs("exome")
#' sigs_denorm <- apply(sigs, 1, function(sig) { 
#'      tmp <- sig[,1] * freqs
#'      tmp / sum(tmp)
#' })
#' @useDynLib sigfit, .registration = TRUE
#' @export
human_trinuc_freqs <- function(type = "genome") {
    if (type == "genome") {
        # Human genome trinucleotide frequencies (from EMu)
        c(1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>A @ AC[ACGT]
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
          1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08) # T>G @ AT[ACGT]
    }
    else if (type == "exome") {
        # Human exome trinucleotide frequencies (from EMu)
        c(1940794, 1442408, 514826, 1403756,
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
          1391660, 1674368, 1559846, 2850934)
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

gen_bar_plot <- function(samples, featurename, title, prob, thresh, 
                         primary_col = "dodgerblue3", secondary_col = "grey90",
                         ...) {
    feature <- rstan::extract(samples, pars = featurename)[[featurename]]
    if ( length(dim(feature)) > 2 && dim(feature)[2] > 1) {
        stop("Plotting for multiple samples not implemented")
    }
    feature <- feature[,1,]
    mean_feature <- colMeans(feature)
    names(mean_feature) <- 1:length(mean_feature)
    error <- HPDinterval(as.mcmc(feature), prob = prob)
    bars <- barplot(mean_feature, ylim = c(0, max(error[, 2]*1.05)), main = title,
                    col = ifelse(error[, 1] > thresh, primary_col, secondary_col), ...)
    top_arr <- arrows(bars, error[, 1], bars, mean_feature, angle=90, code=1, length=0.05)
    bottom_arr <- arrows(bars, error[, 2], bars, mean_feature, angle=90, code=1, length=0.05)
    list(bars, top = top_arr, bottom = bottom_arr)
}

#' Plots estimated exposure of each signature
#' @param samples The MCMC samples
#' @param prob The width of the HPD interval
#' @param thresh Signatures with a lower HPDI below this value are coloured grey
#' @param title The main title of this plot
#' @param primary_col The colour to plot strongly supported exposures (default=dodgerblue3)
#' @param secondary_col The colour to plot weakly supported exposures (default=grey90)
#' @param ... Arguments to pass through to graphics::barplot
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @importFrom "rstan" summary
#' @importFrom "rstan" extract
#' @export
plot_exposures <- function(samples, prob = 0.9, thresh = 1e-3,
                           title = "Signature exposures",
                           primary_col = "dodgerblue3", secondary_col = "grey90",
                           ...) {
    plt <- gen_bar_plot(samples, "exposures", title, prob, thresh, 
                        primary_col, secondary_col, ...)
    plt$bars
    plt$top
    plt$bottom
    legend("topright", legend = sprintf("%.0f%% HPD > %.3f", prob*100, thresh), fill = primary_col)
}

# Plots the fitted spectrum
# @param samples The MCMC samples
# @param prob The width of the HPD interval
# @param title The main title of this plot
# @importFrom "rstan" extract
# @export
# plot_spectrum <- function(samples, prob = 0.9, title = "Fitted spectrum", ...) {
#     plt <- gen_bar_plot(samples, "probs", title, prob, 0, ...)
#     plt$bars
#     plt$top
#     plt$bottom
# }

#' Plots one or more spectra, which can be either mutational catalogues or mutational
#' signatures. For multiple spectra, one plot per spectrum will be generated.
#' @spectra Either a vector with 96 elements, a matrix with 96 columns and one row per signature/catalogue,
#' or a list of signature matrices, obtained via retrieve_pars(). In the latter case, HPD intervals will
#' also be plotted.
#' @counts Do the values in 'spectra' represent mutation counts instead of mutation probabilities?
#' @examples
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "emu", 
#' opportunities = "human-genome")
#' 
#' # Retrieve array of signatures
#' signatures <- retrieve_pars(samples, "signatures")
#' 
#' # Plot mean exposures
#' plot_spectrum(signatures)
#' @useDynLib sigfit, .registration = TRUE
#' @export
plot_spectrum <- function(spectra, counts = FALSE) {
    NCAT <- 96  # number of categories
    # Fetch HPD interval values, if present
    if (is.list(spectra)) {
        spec <- spectra[[1]]
        lwr <- spectra[[2]]
        upr <- spectra[[3]]
    }
    else {
        spec <- spectra
        lwr <- NULL
        upr <- NULL
    }
    # Ensure spectrum is a matrix (96 columns)
    if (is.vector(spec)) 
        spec <- matrix(spec, nrow = 1)
    stopifnot(ncol(spec) == NCAT)
    
    # Plot each spectrum
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    TYPES <- c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
    max_y <- ifelse(is.null(upr), max(spec) * 1.1, max(upr) * 1.1)
    if (!counts)
        probs <- seq(0, 1, 0.05)
        max_y <- probs[which.max(probs > max_y)]
    for (i in 1:nrow(spec)) {
        # Plot spectrum bars
        bars <- barplot(spec[i,],
                        col = rep(COLORS, each = 16), border = "white",
                        yaxt = "n", ylim = c(0, max_y), xlim = c(-1, 116),
                        cex = 1.3, cex.axis = 1.5, cex.lab = 1.7, 
                        las = 2, xaxs = "i", family = "mono")
        axis(side = 2, at = seq(0, max_y, 0.05), las = 2, cex.axis = 1.25)
        label <- ifelse(counts, "Mutation count", "Mutation probability")
        n_text <- ifelse(counts, paste0(" (N = ", sum(spec[i,]), " mutations)"), "")
        mtext(label, side = 2, cex = 1.7, line = 4.5)
        title(paste0("Mutational spectrum #", i, "\n", rownames(spec)[i], n_text), 
              line = 1.5, cex.main = 2)
        # Plot HPD intervals
        if (!is.null(lwr)) {
            arrows(bars, spec[i,], bars, lwr[i,], angle = 90, 
                   length = 0.05, lwd = 1.5, col = "darkgrey")
            arrows(bars, spec[i,], bars, upr[i,], angle = 90, 
                   length = 0.05, lwd = 1.5, col = "darkgrey")
        }
        # Plot mutation type labels
        xleft <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
        xright <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
        rect(xleft = xleft, xright = xright, ybottom = max_y * 0.95, ytop = max_y, 
             col = COLORS, border = "white") #, lty = par("lty"), lwd = par("lwd"))
        text(x = (xleft + xright) / 2, y = max_y * 0.9, labels = TYPES, cex = 2.25)
    }
}

#' Obtains summary values for a set of model parameters (signatures or exposures) from a stanfit object.
#' @param object An object of class stanfit.
#' @param feature Name of the parameter set to extract; either "signatures" or "exposures".
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the HPD intervals.
#' @param signature_names Vector containing the names of the signatures used for fitting. Used only when 
#' retrieving exposures from fitted signatures.
#' @examples
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "emu", 
#' opportunities = "human-genome")
#' 
#' # Retrieve array of signatures
#' signatures <- retrieve_pars(samples, "signatures")
#' 
#' # Retrieve array of exposures using 90% HPD intervals
#' exposures <- retrieve_pars(samples, "exposures", prob = 0.9)
#' 
#' # Plot mean exposures
#' barplot(exposures$mean)  ## or barplot(exposures[[1]])
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @export
retrieve_pars <- function(object, feature, prob = 0.95, signature_names = NULL) {
    feat <- extract(object, pars = feature)[[feature]]
    # Multi-sample case
    if (length(dim(feat)) > 2) {
        # Assign dimension names
        if (feature == "signatures") {
            names2 <- mut_types()
            if (is.null(signature_names)) {
                names1 <- paste("Signature", LETTERS[1:dim(feat)[2]])
            }
            else {
                if (dim(feat)[2] != length(signature_names)) 
                    stop("signature_names must have length equal to number of signatures")
                names1 <- signature_names
            }
        }
        else if (feature == "exposures") {
            names1 <- NULL
            if (is.null(signature_names)) {
                names2 <- paste("Signature", LETTERS[1:dim(feat)[3]])
            }
            else {
                if (dim(feat)[3] != length(signature_names)) 
                    stop("signature_names must have length equal to number of signatures")
                names2 <- signature_names
            }
        }
        # for signatures: Signatures x Categories matrix
        # for exposures: Samples x Signatures matrix
        feat.summ <- list(matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)),
                          matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)),
                          matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)))
        names(feat.summ) <- c("mean", paste0(c("lower_", "upper_"), prob * 100))
        for (i in 1:nrow(feat.summ[[1]])) {
            hpd <- HPDinterval(as.mcmc(feat[,i,]), prob = prob)
            feat.summ[[1]][i,] <- colMeans(feat[,i,])
            feat.summ[[2]][i,] <- hpd[,1]
            feat.summ[[3]][i,] <- hpd[,2]
        }
    } 
    # Single-sample case (only possible in fitting)
    else {
        names1 <- NULL
        if (is.null(signature_names)) {
            names2 <- paste("Signature", LETTERS[1:dim(feat)[3]])
        }
        else {
            if (dim(feat)[3] != length(signature_names)) 
                stop("signature_names must have length equal to number of signatures")
            names2 <- signature_names
        }
        feat.summ <- list(matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)),
                          matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)),
                          matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)))
        names(feat.summ) <- c("mean", paste0(c("lower_", "upper_"), prob * 100))
        for (i in 1:nrow(feat.summ[[1]])) {
            hpd <- HPDinterval(as.mcmc(feat), prob = prob)
            feat.summ[[1]][1,] <- colMeans(feat)
            feat.summ[[2]][1,] <- hpd[,1]
            feat.summ[[3]][1,] <- hpd[,2]
        }
    }
    feat.summ
}

#' Runs MCMC to fit signatures and estimate exposures
#' @param counts Matrix of mutation counts per category (columns) per genome sample (rows).
#' @param signatures Matrix of mutational signatures (rows) to be fitted.
#' @param prior Vector of the same length as signatures, to be used as the Dirichlet prior in the sampling chain. Default prior is uniform (uninformative).
#' @param method Either "nmf" (default) or "emu".
#' @param opportunities Optional matrix of mutational opportunities for "emu" method; must have same dimension as counts. 
#' If equals to "human-genome" or "human-exome", the reference human genome/exome opportunities will be used for every sample.
#' @param ... Arguments to pass to rstan::sampling.
#' @examples
#'  # Custom prior favours signature 1 over 2, 3 and 4
#' samples <- fit_signatures(mycounts, mysignatures, prior = c(5, 1, 1, 1))
#' 
#' # Run a single chain for quite a long time
#' samples <- fit_signatures(mycounts, mysignatures, chains = 1, niter = 13000, warmup = 3000)
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @export
fit_signatures <- function(counts, signatures, prior = NULL,
                           method = "nmf", opportunities = NULL, ...) {
    if (is.null(prior)) {
        prior = rep(1, nrow(signatures))
    }
    
    # If counts is a vector, convert it to a single row matrix
    if (is.vector(counts)) counts <- matrix(counts, nrow = 1)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NSIG, NCAT]
    stopifnot(ncol(counts) == ncol(signatures))
    stopifnot(length(prior) == nrow(signatures))
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    if (method == "emu") {
        if (is.null(opportunities)) {
            opportunities <- matrix(1, nrow = nrow(counts), ncol = ncol(counts))
        }
        else if (opportunities == "human-genome") {
            opportunities <- matrix(rep(human_trinuc_freqs("genome"), nrow(counts)),
                                    nrow = nrow(counts), ncol = ncol(counts), byrow = T)
        }
        else if (opportunities == "human-exome") {
            opportunities <- matrix(rep(human_trinuc_freqs("exome"), nrow(counts)),
                                    nrow = nrow(counts), ncol = ncol(counts), byrow = T)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        # NEED TO IMPLEMENT alpha
        dat <- list(
            C = ncol(counts),
            S = nrow(signatures),
            G = nrow(counts),
            signatures = as.matrix(signatures),
            counts = as.matrix(counts),
            opps = as.matrix(opportunities),
            alpha = rep(1, nrow(signatures))  # TODO: properly build/pass alpha vector
        )
        model <- stanmodels$sigfit_fit_emu
    }
    
    else {
        dat <- list(
            C = ncol(counts),
            S = nrow(signatures),
            G = nrow(counts),
            signatures = t(as.matrix(signatures)),
            counts = as.matrix(counts),
            alpha = prior
        )
        model <- stanmodels$sigfit_fit_nmf
    }
    sampling(model, data = dat, ...)
}

#' Extracts signatures from a set of mutation counts
#' using models based on NMF or EMu
#' 
#' @param counts Matrix of mutation counts for each sample (rows) in each category
#' (columns).
#' @param nsignatures Number (or range of numbers) of signatures to extract.
#' @param method Either "emu" (default) or "nmf" (though currently "nmf" is experimental).
#' @param opportunities Optional matrix of mutational opportunities for "emu" method; must have same dimension as counts. 
#' If equals to "human-genome" or "human-exome", the reference human genome/exome opportunities will be used for every sample.
#' @param exposures_prior Single numeric value to use for all the exposure priors; by default, 0.5 (i.e. Jeffreys prior).
#' @param stanfunc "sampling"|"optimizing"|"vb" Choice of rstan inference strategy. 
#' "sampling" is the full Bayesian MCMC approach, and is the default. "optimizing"
#' returns the Maximum a Posteriori (MAP) point estimates via numerical optimization.
#' "vb" uses Variational Bayes to approximate the full posterior.
#' @param ... Any other parameters to pass through to rstan.
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @importFrom "loo" loo extract_log_lik
#' @export
extract_signatures <- function(counts, nsignatures, method = "emu", 
                               opportunities = NULL, exposures_prior = 0.5, 
                               stanfunc = "sampling", ...) {
    
    if (method == "emu") {
        if (is.null(opportunities)) {
            opportunities <- matrix(1, nrow = nrow(counts), ncol = ncol(counts))
        }
        else if (opportunities == "human-genome") {
            opportunities <- matrix(rep(human_trinuc_freqs("genome"), nrow(counts)),
                                    nrow = nrow(counts), ncol = ncol(counts), byrow = TRUE)
        }
        else if (opportunities == "human-exome") {
            opportunities <- matrix(rep(human_trinuc_freqs("exome"), nrow(counts)),
                                    nrow = nrow(counts), ncol = ncol(counts), byrow = TRUE)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        model <- stanmodels$sigfit_ext_emu
        dat <- list(
            C = ncol(counts),
            G = nrow(counts),
            S = nsignatures[1],
            counts = as.matrix(counts),
            opps = as.matrix(opportunities)
        )
    }
    else if (method == "nmf") {
        if (!is.null(opportunities)) {
            warning("Using \"nmf\" model: opportunities will be ignored.")
        }
        
        model <- stanmodels$sigfit_ext_nmf
        dat <- list(
            C = ncol(counts),
            G = nrow(counts),
            S = nsignatures[1],
            counts = as.matrix(counts),
            exposures_prior_val = exposures_prior
        )
    }
    
    # Extract signatures for each nsignatures value
    if (length(nsignatures) > 1) {
        out <- vector(mode = "list", length = max(nsignatures))
        for (n in nsignatures) {
            cat("Extracting", n, "signatures\n")
            dat$S <- as.integer(n)
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
        
        # Identify best number of signatures using LOOIC
        out$looic <- rep(NA, max(nsignatures))
        for (n in nsignatures) {
            out$looic[n] <- loo(loo::extract_log_lik(out[[n]]))$looic
        }
        out$best_N <- which.min(out$looic)
        cat("Best number of signatures (lowest LOOIC) is", out$best_N)
    }
    # Single nsignatures value case
    else {
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
    # counts[NSAMPLES, NCAT], signatures[NSIG, NCAT]
    stopifnot(ncol(counts) == ncol(signatures))
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    dat <- list(
        C = ncol(counts),
        S = nrow(signatures),
        G = nrow(counts),
        N = num_extra_sigs,
        fixed_sigs = as.matrix(signatures),
        counts = as.matrix(counts)
    )
    ## Only NMF implemented so far
    model <- stanmodels$sigfit_fitex_nmf
    
    if (stanfunc == "sampling") {
        cat("Stan sampling:")
        sampling(model, data = dat, chains = 1, ...)
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
