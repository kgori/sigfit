#' Find the best matches between two sets of signatures
#' 
#' \code{match_signatures} compares two independent estimates of signatures to 
#' find the closest matches between them.
#' @param sigs_a Signatures estimate: either a list such as one produced by 
#' \code{\link{retrieve_pars}}, with a \code{$mean} entry, or a matrix with one 
#' row per signature and one column for each of the 96 mutation types.
#' @param sigs_b Signatures estimate as for \code{sigs_a}.
#' @return A vector containing, for each signature in \code{sigs_a}, the index 
#' of the closest match in \code{sigs_b}.
#' @importFrom "clue" solve_LSAP
#' @export
match_signatures <- function(sigs_a, sigs_b) {
    if ("mean" %in% names(sigs_a)) a <- sigs_a$mean
    else a <- sigs_a
    
    if ("mean" %in% names(sigs_b)) b <- sigs_b$mean
    else b <- sigs_b
    
    nA <- nrow(a)
    nB <- nrow(b)
    stopifnot(nA == nB)
    
    m <- matrix(0.0, nrow = nA, ncol = nB)
    for (i in 1:nA) {
        for (j in 1:nB) {
            m[i, j] <- cosine_sim(a[i, ], b[j, ])
        }
    }
    solve_LSAP(m, maximum = TRUE)
}

#' Initial coercion to matrix for signatures/exposures/counts
to_matrix <- function(x) {
    # If x is coming from retrieve_pars, get mean
    if (is.list(x) & "mean" %in% names(x))  x <- x$mean
    # If x is a vector, transform to 1-row matrix
    if (is.vector(x))  x <- matrix(x, nrow = 1)
    # Otherwise, try coercing to matrix
    if (!is.matrix(x))  x <- as.matrix(x)
    x
}

#' Build a opportunities matrix
build_opps_matrix <- function(nsamples, keyword, strand) {
    matrix(rep(human_trinuc_freqs(keyword, strand), nsamples),
           nrow = nsamples, 
           byrow = TRUE)
}

#' Cosine similarity between two vectors
cosine_sim <- function(x, y) { x %*% y / sqrt(x%*%x * y%*%y) }

#' L2 norm between two vectors
l2_norm <- function(x, y) { sqrt(sum((x - y)^2)) }

#' Deal with zero values in a signature
remove_zeros_ <- function(mtx, min_allowed = 1e-9) {
    as.matrix(
        t(apply(mtx, 1, function(row) {
            row[row < min_allowed] <- min_allowed
            row / sum(row)
        }))
    )
}

#' From https://stackoverflow.com/a/21689613/1875814
#' Produce Excel-style letter labels - A, B, ..., Z, AA, AB, ..., AZ, ..., ZZ, AAA, ...
letterwrap <- function(n, depth = 1) {
    args <- lapply(1:depth, FUN = function(x) return(LETTERS))
    x <- do.call(expand.grid, args = list(args, stringsAsFactors = F))
    x <- x[, rev(names(x)), drop = F]
    x <- do.call(paste0, x)
    if (n <= length(x)) return(x[1:n])
    return(c(x, letterwrap(n - length(x), depth = depth + 1)))
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

#' Generate default colour palette for signatures
default_sig_palette <- function(n) {
    rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
          "#A6761D", "#666666", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
          "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"), 
        5)[1:n]
}

#' Generates character vector of 96 (or 192, if strand=T) trinucleotide mutation types
mut_types <- function(strand = FALSE) {
    bases <- c("A", "C", "G", "T")
    muts <- paste0(rep(rep(bases, each = 4), 6),
                   rep(bases[c(2, 4)], each = 48),
                   rep(bases, 6 * 16 / 4),
                   ">",
                   rep(rep(bases, each = 4), 6),
                   c(rep(bases[-2], each = 16), rep(bases[-4], each = 16)),
                   rep(bases, 6 * 16 / 4))
    if (strand) {
        paste(c(rep("T", 96), rep("U", 96)), muts, sep = ":")
    }
    else {
        muts
    }
}

#' Access stan models
#' @export
stan_models <- function() {
    stanmodels
}

#' Retrieve human trinucleotide frequencies
#' 
#' \code{human_trinuc_freqs} returns the reference human genome or exome 
#' trinucleotide frequencies.
#' @param type Character; either \code{"genome"} (default) or \code{"exome"}.
#' @param strand Logical; if \code{TRUE}, a strand-bias representation of catalogues
#' and signatures will be used.
#' @return A numeric vector containinig 96 frequency values (one per trinucleotide type), if
#' \code{strand=FALSE}, or 192 frequency values (one per trinucleotide and strand type), if
#' \code{strand=TRUE}. In the latter case, the trinucleotide frequencies are assumed to be
#' equally distributed between the two strands.
#' @export
human_trinuc_freqs <- function(type = "genome", strand = FALSE) {
    if (type == "genome") {
        # Human genome trinucleotide frequencies (from EMu)
        freq <- c(1.14e+08, 6.60e+07, 1.43e+07, 9.12e+07, # C>A @ AC[ACGT]
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
        freq <- c(1940794, 1442408, 514826, 1403756,
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
        stop("'type' must be either \"genome\" or \"exome\"")
    }
    
    if (strand) {
        rep(freq / 2, 2)
    }
    else {
        freq
    }
}

#' Build mutational catalogues
#' 
#' \code{build_catalogues} generates a set of mutational catalogues from a table containing 
#' the base change and trinucleotide context of each single-nucleotide variant in every sample.
#' @param variants Matrix with one row per single-nucleotide variant, and four or five columns:
#' \itemize{
#'  \item{Sample ID (character, e.g. "Sample 1").}
#'  \item{Reference allele (character: "A", "C", "G", or "T").}
#'  \item{Alternate allele (character: "A", "C", "G", or "T").}
#'  \item{Trinucleotide context of the variant (the reference sequence between the positions 
#'  immediately before and after the variant; character, e.g. "TCA").}
#'  \item{Optional: transcriptional strand of the variant (character: "T" for transcribed,
#'  or "U" for untranscribed). If this column is included, a strand-bias representation of
#'  catalogues will be used.}
#' }
#' @return A matrix of mutation counts, where each row corresponds to a sample and each column
#' corresponds to one of the 96 trinucleotide mutation types.
#' @examples
#' # Load example mutation data
#' data("variants_21breast")
#' head(variants_21breast)
#' 
#' # Build catalogues
#' counts <- build_catalogues(variants_21breast)
#' counts
#' @export
build_catalogues <- function(variants) {
    if (!(ncol(variants) %in% c(4, 5))) {
        stop("'variants' must be a matrix with 4 or 5 columns and one row per variant.\nType ?build_catalogues to see the documentation.")
    }
    if (ncol(variants) == 5) {
        if (!all(unique(variants[, 5] %in% c("T", "U")))) {
            stop("The fifth column of 'variants' (transcriptional strand) can only contain \"T\" or \"U\" values.\nType ?build_catalogues to see the documentation.")
        }
        cat("'variants' has 5 columns: strand-bias catalogues will be generated\n")
        strand <- TRUE
    }
    else {
        strand <- FALSE
    }
    
    # Check that REF base coincides with middle base in trinucleotide
    if (any(variants[, 2] != substr(variants[, 4], 2, 2))) {
        stop("REF base (column 2) is not equal to middle base of the trinucleotide (column 4).")
    }
    
    # Make catalogue matrix with one row per sample
    samples <- unique(variants[, 1])
    catalogues <- t(sapply(samples, function(sample) {
        # Select mutations from sample
        idx <- variants[, 1] == sample
        
        # Obtain mutation types, collapsed such that they refer to pyrimidine bases
        vars_collapsed <- apply(variants[idx, ], 1, function(var) {
            if (var[2] %in% c("C", "T")) {
                trinuc <- strsplit(var[4], "")[[1]]
                alt <- var[3]
                if (strand) {
                    strd <- var[5]
                }
            }
            else {
                trinuc <- rev_comp(var[4])
                alt <- rev_comp(var[3])
                if (strand) {
                    strd <- ifelse(var[5] == "U", "T", "U")
                }
            }
            if (strand) {
                paste0(strd, ":", paste(trinuc, collapse=""), ">", trinuc[1], alt, trinuc[3])
            }
            else {
                paste0(paste(trinuc, collapse=""), ">", trinuc[1], alt, trinuc[3])
            }
        })
        
        # Count number of occurrences of each mutation type
        stopifnot(all(vars_collapsed %in% mut_types(strand)))
        sapply(mut_types(strand), function(type) {
            sum(grepl(type, vars_collapsed, fixed = TRUE))
        })
    }))
}

#' Fetch COSMIC mutational signatures
#' 
#' \code{fetch_cosmic_data} downloads the latest release of signatures from COSMIC
#' (http://cancer.sanger.ac.uk/cosmic/signatures) and produces a matrix of signatures
#' that can be used for signature fitting.
#' @param reorder If \code{TRUE}, the matrix will be reordered by substitution type and trinucleotide.
#' @param remove_zeros If \code{TRUE}, pseudocounts will be added to prevent the signatures from
#' containing any zeros, which can affect computation of the log likelihood.
#' @return Matrix of signatures, with one row per signature and one column for each of
#' the 96 trinucleotide mutation types.
#' @importFrom "utils" read.table
#' @export
fetch_cosmic_data <- function(reorder = TRUE, remove_zeros = TRUE) {
    cosmic_sigs <- read.table("http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", 
                              header = TRUE, sep = "\t", check.names = FALSE)
    if (reorder) {
        cosmic_sigs <- cosmic_sigs[order(cosmic_sigs[["Substitution Type"]], cosmic_sigs[["Trinucleotide"]]),]
    }
    rownames(cosmic_sigs) <- mut_types()
    cosmic_sigs <- t(cosmic_sigs[, paste("Signature", 1:30)])
    if (remove_zeros) {
        cosmic_sigs <- remove_zeros_(cosmic_sigs)
    }
    cosmic_sigs
}

#' Convert signatures between models
#' 
#' \code{convert_signatures} converts between the representation of signatures used
#' in the NMF model (which is relative to the reference mutational opportunities), and
#' the representation used in the EMu model (which is not relative to mutational opportunities).
#' This is done by multiplying or dividing each signature by the average mutational opportunities
#' of the samples, or by the human genome/exome reference trinucleotide frequencies.
#' @param signatures Either a matrix of mutational signatures, with one row per signature and one
#' column for each of the 96 (or 192) mutation types, or a list of signatures generated via
#' \code{\link{retrieve_pars}}.
#' @param ref_opportunities Numeric vector of reference or average mutational opportunities, with 
#' one element for each of the 96 (or 192) mutation types. If equal to \code{"human-genome"} or 
#' \code{"human-exome"}, the reference human genome/exome mutational opportunities will be used.
#' @param model_to The model to convert to: either \code{"nmf"} (in which case the signatures will
#' be multiplied by the opportunities) or \code{"emu"} (in which case the signatures will be divided
#' by the opportunities).
#' @return A matrix of transformed signatures with the same dimensions as \code{signatures}.
#' @examples
#' # Fetch COSMIC signatures 
#' # These are in "NMF" format, i.e. they are relative
#' # to the human genome mutational opportunities
#' signatures <- fetch_cosmic_data()
#' 
#' # Plot COSMIC signature 1
#' barplot(signatures[1,])
#' 
#' # Convert signatures to the "EMu" format, i.e. make 
#' # them not relative to mutational opportunities
#' converted_signatures <- convert_signatures(signatures,
#'                                            ref_opportunities = "human-genome",
#'                                            model_to = "emu")
#' 
#' # Plot COSMIC signature 1, converted to "EMu" format
#' barplot(converted_signatures[1,])
#' @export
convert_signatures <- function(signatures, ref_opportunities, model_to) {
    signatures <- to_matrix(signatures)
    stopifnot(ncol(signatures) %in% c(96, 192))
    strand <- ncol(signatures) == 192
    if (strand) {
        cat("'signatures' contains 192 mutations types: using strand-bias opportunities.\n")
    }
    
    if (ref_opportunities == "human-genome") {
        ref_opportunities <- human_trinuc_freqs(type = "genome", strand = strand)
    }
    else if (ref_opportunities == "human-exome") {
        ref_opportunities <- human_trinuc_freqs(type = "exome", strand = strand)
    }
    ref_opportunities <- as.numeric(ref_opportunities)
    if (length(ref_opportunities) != ncol(signatures)) {
        stop("'ref_opportunities' must have one element for each column in 'signatures'.")
    }
    
    if (model_to == "nmf") {
        t(apply(signatures, 1, function(row) {
            x <- row * ref_opportunities
            x / sum(x)
        }))
    }
    else if (model_to == "emu") {
        t(apply(signatures, 1, function(row) {
            x <- row / ref_opportunities
            x / sum(x)
        }))
    }
    else {
        stop("'model_to' must be either \"nmf\" or \"emu\".")
    }
}

#' Retrieve model parameters
#' 
#' Obtains summary values for a set of model parameters (signatures or exposures) from a stanfit object.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}.
#' @param feature Name of the parameter set to extract; either \code{"signatures"} or \code{"exposures"}.
#' @param hpd_prob A value in the interval (0, 1), giving the target probability content of 
#' the HPD intervals.
#' @param signature_names Vector containing the names of the signatures used for fitting. Used only when 
#' retrieving exposures from fitted signatures.
#' @return A list of three matrices, which respectively contain the values corresponding to the
#' mean of the model parameter of interest, and to the lower and upper ends of its HPD interval.
#' @examples
#' # Load example mutational catalogues
#' data("counts_21breast")
#' 
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(counts_21breast, nsignatures = 2, method = "emu",
#'                               opportunities = "human-genome", iter = 800)
#' 
#' # Retrieve signatures and exposures
#' signatures <- retrieve_pars(samples, "signatures")
#' 
#' # Retrieve array of exposures using 90% HPD intervals
#' exposures <- retrieve_pars(samples, "exposures", hpd_prob = 0.9)
#' 
#' # Plot signatures
#' plot_spectrum(signatures)
#' 
#' # Plot mean exposures
#' barplot(exposures$mean)   ## or: barplot(exposures[[1]])
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @export
retrieve_pars <- function(mcmc_samples, feature, hpd_prob = 0.95, counts = NULL, signatures = NULL, signature_names = NULL) {
    if (feature == "reconstructions") {
        if (is.null(counts)) {
            stop("Counts must be provided to retrieve reconstructions")
        }
        l <- get_reconstructions(counts, mcmc_samples, signatures)
        feat_summ <- list(mean = apply(l$reconstructions, c(1, 3), sum),
                          lower = l$hpds[, 1, ],
                          upper = l$hpds[, 2, ])
    }
    else {
        feat <- extract(mcmc_samples, pars = feature)[[feature]]
        strand <- dim(feat)[3] == 192  # strand bias indicator
        
        # Multi-sample case
        if (length(dim(feat)) > 2) {
            # Assign dimension names
            if (feature == "signatures") {
                names2 <- mut_types(strand)
                if (is.null(signature_names)) {
                    LETTERLABELS <- letterwrap(dim(feat)[2])
                    names1 <- paste("Signature", LETTERLABELS[1:dim(feat)[2]])
                }
                else {
                    if (dim(feat)[2] != length(signature_names)) {
                        stop("'signature_names' must have length equal to the number of signatures.")
                    }
                    names1 <- signature_names
                }
            }
            else if (feature == "exposures") {
                names1 <- NULL
                if (is.null(signature_names)) {
                    LETTERLABELS <- letterwrap(dim(feat)[3])
                    names2 <- paste("Signature", LETTERLABELS[1:dim(feat)[3]])
                }
                else {
                    if (dim(feat)[3] != length(signature_names))  {
                        stop("'signature_names' must have length equal to the number of signatures")
                    }
                    names2 <- signature_names
                }
            }
            # for signatures: Signatures x Categories matrix
            # for exposures: Samples x Signatures matrix
            feat_summ <- list(matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)),
                              matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)),
                              matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)))
            names(feat_summ) <- c("mean", paste0(c("lower_", "upper_"), hpd_prob * 100))
            for (i in 1:nrow(feat_summ[[1]])) {
                hpd <- HPDinterval(as.mcmc(feat[,i,]), prob = hpd_prob)
                feat_summ[[1]][i,] <- colMeans(feat[,i,])
                feat_summ[[2]][i,] <- hpd[,1]
                feat_summ[[3]][i,] <- hpd[,2]
            }
        }
        
        # Single-sample case (only possible when fitting)
        else {
            names1 <- NULL
            if (is.null(signature_names)) {
                LETTERLABELS <- letterwrap(dim(feat)[3])
                names2 <- paste("Signature", LETTERLABELS[1:dim(feat)[3]])
            }
            else {
                if (dim(feat)[3] != length(signature_names)) {
                    stop("'signature_names' must have length equal to the number of signatures.")
                }
                names2 <- signature_names
            }
            feat_summ <- list(matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)),
                              matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)),
                              matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)))
            names(feat_summ) <- c("mean", paste0(c("lower_", "upper_"), hpd_prob * 100))
            for (i in 1:nrow(feat_summ[[1]])) {
                hpd <- HPDinterval(as.mcmc(feat), prob = hpd_prob)
                feat_summ[[1]][1,] <- colMeans(feat)
                feat_summ[[2]][1,] <- hpd[,1]
                feat_summ[[3]][1,] <- hpd[,2]
            }
        }
    }
    feat_summ
}

#' Generate posterior predictive check values from a model
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}.
#' @importFrom "stats" rmultinom rpois
#' @export
simulate_ppc <- function(counts, mcmc_samples) {
    if (grepl("nmf", mcmc_samples@model_name)) {
        e <- extract(mcmc_samples, pars = "probs")
        ppc <- array(0, dim(e$probs))
        for (i in 1:dim(e$probs)[1]) {
            for (j in 1:dim(e$probs)[2]) {
                ppc[i, j, ] <- t(
                    rmultinom(1, sum(counts[j, ]), e$probs[i, j, ])
                )
            }
        }
    }
    else {
        e <- extract(mcmc_samples, pars = "expected_counts")
        ppc <- array(0, dim(e$expected_counts))
        for (i in 1:dim(e$expected_counts)[1]) {
            for (j in 1:dim(e$expected_counts)[2]) {
                for (k in 1:dim(e$expected_counts)[3]) {
                    ppc[i, j, k] <- rpois(1, e$expected_counts[i, j, k])
                }
            }
        }
    }
    ppc
}

#' Generate log likelihood values from a model
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}.
#' @importFrom "stats" dmultinom dpois
#' @export
get_loglik <- function(counts, mcmc_samples) {
    dnames <- list(NULL, NULL)
    names(dnames) <- c("iterations", "")
    if (grepl("nmf", mcmc_samples@model_name)) {
        e <- extract(mcmc_samples, pars = "probs")
        nrep <- dim(e$probs)[1]
        nsamples <- dim(e$probs)[2]
        log_lik <- matrix(0, nrow = nrep, ncol = nsamples, dimnames = dnames)
        for (i in 1:nrep) {
            for (j in 1:nsamples) {
                log_lik[i, j] <- dmultinom(counts[j, ], prob = e$probs[i, j, ], log = TRUE)
            }
        }
    }
    else {
        e <- extract(mcmc_samples, pars = "expected_counts")
        nrep <- dim(e$expected_counts)[1]
        nsamples <- dim(e$expected_counts)[2]
        log_lik <- matrix(0, nrow = nrep, ncol = nsamples, dimnames = dnames)
        for (i in 1:nrep) {
            for (j in 1:nsamples) {
                log_lik[i, j] <- sum(dpois(counts[j, ], e$expected_counts[i, j, ], log = TRUE))
            }
        }
    }
    log_lik
}

#' Generate reconstructed mutation catalogues from parameters estimated from MCMC samples
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}.
#' @param signatures Matrix of signatures to use for fitting models, which make no estimate of signatures.
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
#' @export
get_reconstructions <- function(counts, mcmc_samples, signatures = NULL) {
    NCAT <- ncol(counts)   # number of categories
    NSAMP <- nrow(counts)  # number of samples
    
    e <- extract(mcmc_samples)
    NREP <- dim(e$exposures)[1]
    stopifnot(NSAMP == dim(e$exposures)[2])
    NSIG <- dim(e$exposures)[3]

    if (!"signatures" %in% names(e)) {
        if (is.null(signatures)) {
            stop("No signatures found in MCMC samples and no matrix of signatures provided")
        }
        e$signatures <- aperm(
            sapply(1:NREP, function(i) as.matrix(signatures), simplify = "array"),
            c(3, 1, 2)
        )
    }
    
    reconstructions <- array(NA, dim = c(NSAMP, NSIG, NCAT))
    hpds <- array(NA, dim = c(NSAMP, 2, NCAT))
    for (sample in 1:NSAMP) {
        if (grepl("emu", mcmc_samples@model_name)) {        
            arr <- aperm(     
                sapply(1:NREP, function(i) {      
                    sweep(e$activities[i, sample, ] * e$signatures[i, , ],
                          2, as.numeric(opportunities[sample, ]), `*`)      
                }, simplify = "array"),       
                c(3, 1, 2)        
            )     
        }     
        
        # For NMF results     
        else {        
            arr <- aperm(     
                sapply(1:NREP, function(i) {      
                    e$exposures[i, sample, ] *        
                        e$signatures[i, , ] *         
                        sum(counts[sample, ])     
                }, simplify = "array"),       
                c(3, 1, 2)        
            )     
        }
        reconstructions[sample, , ] <- apply(arr, c(2, 3), mean)
        hpds[sample, , ] <- t(HPDinterval(
            as.mcmc(apply(arr, c(1, 3), sum))
        ))
    }
    list(reconstructions = reconstructions, hpds = hpds, exposures = t(apply(e$exposures, 2, colMeans)))
}
