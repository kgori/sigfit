#' Initial coercion to matrix for signatures/exposures/counts
to_matrix <- function(x, int = FALSE) {
    # If x is coming from retrieve_pars, get mean
    if (is.list(x) & "mean" %in% names(x))
        x <- x$mean
    # If x is a vector, transform to 1-row matrix
    if (is.vector(x))
        x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
    # Otherwise, try coercing to matrix
    if (!is.matrix(x))
        x <- as.matrix(x)
    # For counts matrix: if real-valued, round
    if (int) {
        x <- round(x)
    }
    x
}

#' Build an opportunities matrix
build_opps_matrix <- function(nsamples, ncat, opps) {
    if (is.null(opps)) {
        # Default: uniform opportunities
        matrix(1, nrow = nsamples, ncol = ncat)
    }
    else {
        if (is.character(opps)) {
            strand <- ncat == 192
            opps <- matrix(rep(human_trinuc_freqs(opps, strand), nsamples),
                           nrow = nsamples, byrow = TRUE)
        }
        else {
            opps <- as.matrix(opps)
        }
        # Normalise to sum to nsamples (to avoid very large/small values)
        opps / sum(opps) * nsamples
    }
}

#' Cosine similarity between two vectors
#' @param x Vector of real values
#' @param y Vector of real values, same length as \code{x}
#' @return Similarity measure between \code{x} and \code{y}, where 0
#' means no similarity, and 1 means maximum similarity.
#' @examples
#' cosine_sim(c(1,2,3,4), c(4,3,2,1))
cosine_sim <- function(x, y) { x %*% y / sqrt(x%*%x * y%*%y) }

#' L2 norm between two vectors
#' @param x Vector of real values
#' @param y Vector of real values, same length as \code{x}
#' @return Distance measure between \code{x} and \code{y}, where 0
#' means no difference.
#' @examples
#' l2_norm(c(1,2,3,4), c(4,3,2,1))
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
        10)[1:n]
}

#' Generates character vector of 96 (or 192, if strand=TRUE) trinucleotide mutation types
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

#' Generate log likelihood values from a model
#' @param mcmc_samples List with elements \code{`data`} and \code{`results`}, produced via either
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}} or \code{\link{fit_extract_signatures}}.
#' @importFrom "stats" dmultinom dpois
get_loglik <- function(mcmc_samples) {
    counts <- mcmc_samples$data$counts_real
    dnames <- list(NULL, NULL)
    names(dnames) <- c("iterations", "")
    if (grepl("nmf", mcmc_samples$result@model_name)) {
        e <- extract(mcmc_samples$result, pars = "probs")
        NREP <- dim(e$probs)[1]
        NSAMP <- dim(e$probs)[2]
        log_lik <- matrix(0, nrow = NREP, ncol = NSAMP, dimnames = dnames)
        for (i in 1:NREP) {
            for (j in 1:NSAMP) {
                log_lik[i, j] <- dmultinom(counts[j, ], prob = e$probs[i, j, ], log = TRUE)
            }
        }
    }
    else {
        e <- extract(mcmc_samples$result, pars = "expected_counts")
        NREP <- dim(e$expected_counts)[1]
        NSAMP <- dim(e$expected_counts)[2]
        log_lik <- matrix(0, nrow = NREP, ncol = NSAMP, dimnames = dnames)
        for (i in 1:NREP) {
            for (j in 1:NSAMP) {
                log_lik[i, j] <- sum(dpois(counts[j, ], e$expected_counts[i, j, ], log = TRUE))
            }
        }
    }
    log_lik
}

#' Generate reconstructed mutation catalogues from parameters estimated from MCMC samples
#' @param mcmc_samples List with elements \code{`data`} and \code{`results`}, produced via either
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}} or \code{\link{fit_extract_signatures}}.
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
get_reconstructions <- function(mcmc_samples) {
    counts <- mcmc_samples$data$counts_real
    NCAT <- ncol(counts)   # number of categories
    NSAMP <- nrow(counts)  # number of samples
    
    e <- extract(mcmc_samples$result)
    NREP <- dim(e$exposures)[1]
    stopifnot(NSAMP == dim(e$exposures)[2])
    NSIG <- dim(e$exposures)[3]
    
    if (!"signatures" %in% names(e)) {
        signatures <- mcmc_samples$data$signatures
        e$signatures <- aperm(
            sapply(1:NREP, function(i) as.matrix(signatures), simplify = "array"),
            c(3, 1, 2)
        )
    }
    
    reconstructions <- array(NA, dim = c(NSAMP, NSIG, NCAT))
    hpds <- array(NA, dim = c(NSAMP, 2, NCAT))
    opportunities <- mcmc_samples$data$opportunities
    for (sample in 1:NSAMP) {
        if (mcmc_samples$data$family != 1) {
            #if (grepl("emu", mcmc_samples$result@model_name)) {
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
                    probs <- sweep(
                        e$exposures[i, sample, ] *
                            e$signatures[i, , ],
                        2, as.numeric(opportunities[sample, ]), `*`)
                    probs / sum(probs) * sum(counts[sample, ])
                }, simplify = "array"),
                c(3, 1, 2)
            )
        }
        reconstructions[sample, , ] <- apply(arr, c(2, 3), mean)
        
        hpds[sample, , ] <- t(HPDinterval(
            as.mcmc(apply(arr, c(1, 3), sum))
        ))
    }
    if ("signatures" %in% names(mcmc_samples$data)) {
        dimnames(reconstructions)[[2]] <- rownames(mcmc_samples$data$signatures)
    }
    list(reconstructions = reconstructions, hpds = hpds, exposures = t(apply(e$exposures, 2, colMeans)))
}

#' Fetch COSMIC mutational signatures (deprecated)
#'
#' \code{fetch_cosmic_data} downloads the latest release of signatures from COSMIC
#' (http://cancer.sanger.ac.uk/cosmic/signatures) and produces a matrix of signatures
#' that can be used for signature fitting.
#' NB. COSMIC signatures are also available in sigfit via the functions
#' \code{data("cosmic_signatures_v2")}, \code{data("cosmic_signatures_v3") and
#' \code{data("cosmic_signatures_v3_strand")}.
#' @param reorder Logical; if \code{TRUE} (default), the matrix will be reordered by substitution
#' type and trinucleotide.
#' @param remove_zeros Logical; if \code{TRUE} (default), pseudocounts will be added to prevent the
#' signatures from containing any zeros, which can affect computation of the log likelihood.
#' @return Matrix of signatures, with one row per signature and one column for each of
#' the 96 trinucleotide mutation types.
#' @importFrom "utils" read.table
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

#' Access stan models
#' @export
stan_models <- function() {
    stanmodels
}

#' Find the best matches between two sets of signatures
#'
#' \code{match_signatures} compares two independent estimates of signatures to
#' find the closest matches between them.
#' @param sigs_a First set of signatures; either a numeric matrix of signatures,
#' with one row per signature and one column per mutation type, or a list of matrices
#' generated via \code{\link{retrieve_pars}}.
#' @param sigs_b Second set of signatures, as for \code{sigs_a}.
#' @param stat Similarity metric to use when comparing signatures. Admits values \code{"cosine"}
#' (default, cosine similarity) or \code{"L2"} (L2 norm or Euclidean distance).
#' @return A numeric vector containing, for each signature in \code{sigs_a}, the index
#' of the closest match in \code{sigs_b}.
#' @importFrom "clue" solve_LSAP
#' @export
match_signatures <- function(sigs_a, sigs_b, stat = "cosine") {
    a <- to_matrix(sigs_a)
    b <- to_matrix(sigs_b)
    nA <- nrow(a)
    nB <- nrow(b)
    # stopifnot(nA == nB)
    sim_fn <- switch(stat,
                     "cosine" = cosine_sim,
                     "L2" = l2_norm)
    if (is.null(sim_fn)) {
        stop(paste0("'stat' only admits values \"cosine\" and \"L2\".\n",
                    "Type ?match_signatures to read the documentation."))
    }
    
    m <- matrix(0.0, nrow = nA, ncol = nB)
    for (i in 1:nA) {
        for (j in 1:nB) {
            m[i, j] <- sim_fn(a[i, ], b[j, ])
        }
    }
    solve_LSAP(m, maximum = TRUE)
}

#' Retrieve human trinucleotide frequencies
#'
#' \code{human_trinuc_freqs} returns the reference human genome or exome
#' trinucleotide frequencies.
#' @param type Character; admits values \code{"genome"} (default) or \code{"exome"}.
#' @param strand Logical; if \code{TRUE}, transcriptional strand-wise representations of
#' catalogues and signatures will be used (default is \code{FALSE}).
#' @return A numeric vector containing 96 frequency values (one per trinucleotide mutation type),
#' if \code{strand=FALSE}, or 192 frequency values (one per mutation and strand type), if
#' \code{strand=TRUE}. In the latter case, the trinucleotide frequencies are assumed to be
#' equally distributed between the two strands.
#' @export
human_trinuc_freqs <- function(type = "human-genome", strand = FALSE) {
    if (type == "human-genome") {
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
    else if (type == "human-exome") {
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
        stop("'type' must be either \"human-genome\" or \"human-exome\"")
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
#' @param variants Character matrix or data frame, with one row per single-nucleotide variant
#' and four/five columns:
#' \itemize{
#'  \item{Sample ID (character, e.g. "Sample 1").}
#'  \item{Reference allele (character: "A", "C", "G", or "T").}
#'  \item{Alternate allele (character: "A", "C", "G", or "T").}
#'  \item{Trinucleotide context of the variant (character; the reference sequence between the
#'  positions immediately before and after the variant; e.g. "TCA").}
#'  \item{Optional: transcriptional strand of the variant (character/numeric: 1 or "1" or "U" for
#'  untranscribed; -1 or "-1" or "T" for transcribed). If this column is included, a
#'  transcriptional strand-wise representation of catalogues will be used.}
#' }
#' @return An integer matrix of mutation counts, where each row corresponds to a sample and each
#' column corresponds to one of the 96 (or 192) mutation types.
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
        stop("'variants' must be a matrix or data frame with four/five columns and one row per variant.\nType ?build_catalogues to read the documentation.")
    }
    if (ncol(variants) == 5) {
        warning("'variants' has five columns: generating transcriptional-strand-wise catalogues.")
        variants[, 5] <- gsub("1", "U",
                              gsub("-1", "T", as.character(variants[, 5])))
        idx <- !(variants[, 5] %in% c("T", "U"))
        if (any(idx)) {
            warning(sum(idx), " variants have an invalid strand value and have been omitted.\n",
                    "Transcriptional strand (column 5) admits only values 1, -1, \"1\", \"-1\", \"T\" and \"U\".\nType ?build_catalogues to read the documentation.")
            variants <- variants[!idx, ]
        }
        strand <- TRUE
    }
    else {
        strand <- FALSE
    }
    if (!is.matrix(variants)) {
        variants <- as.matrix(variants)
    }

    # Exclude any trinucleotides containing undefined bases
    idx <- !grepl("(A|C|G|T){3}", variants[, 4])
    if (any(idx)) {
        warning(sum(idx), " variants have an invalid trinucleotide context and have been omitted.")
        variants <- variants[!idx, ]
    }

    # Check that REF base coincides with middle base in trinucleotide
    if (any(variants[, 2] != substr(variants[, 4], 2, 2))) {
        stop("Reference allele (column 2) does not match the middle base of the trinucleotide context (column 4) in some variant(s).")
    }

    # Make catalogue matrix with one row per sample
    samples <- unique(variants[, 1])
    catalogues <- t(sapply(samples, function(sample) {
        # Select mutations from sample
        idx <- variants[, 1] == sample

        # Obtain mutation types, collapsed such that they refer to pyrimidine bases
        vars_collapsed <- apply(variants[idx, , drop = FALSE], 1, function(var) {
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
    catalogues
}

#' Convert signatures between models
#'
#' \code{convert_signatures} converts mutational signatures between two representations relative to 
#' different sets of mutational opportunities.
#' @param signatures Either a numeric matrix of mutational signatures, with one row per signature
#' and one column per mutation type, or a list of matrices generated via \code{\link{retrieve_pars}}.
#' @param opportunities_from Mutational opportunities that are currently imposed on the signatures.
#' This can be a numeric matrix of mutational opportunities, with one row per signature and
#' one column per mutation type (if the signatures were inferred together, then every row should
#' contain the same opportunities). Character values \code{"human-genome"} or \code{"human-exome"}
#' are also admitted, in which case the mutational opportunities of the reference human
#' genome/exome will be used for every signature. If this argument is provided, the signatures will
#' be normalised (divided) by the opportunities. Otherwise (if \code{opportunities_from = NULL}),
#' the signatures will be interpreted as being already normalised (i.e. not relative to any
#' opportunities).
#' @param opportunities_to Mutational opportunities that are to be imposed on the signatures.
#' Admits the same values as \code{opportunities_from}. If this argument is provided, the
#' signatures will be multiplied by the opportunities. Otherwise (if \code{opportunities_to = NULL}),
#' the signatures will be converted to a normalised representation (i.e. not relative to any
#' opportunities).
#' @return A numeric matrix of transformed signatures with the same dimensions as \code{signatures}.
#' @examples
#' # Load COSMIC signatures
#' # (these are relative to the human genome mutational opportunities)
#' data("cosmic_signatures_v2")
#'
#' # Plot COSMIC signature 1
#' barplot(cosmic_signatures_v2[1,])
#'
#' # Convert signatures from genome-relative to exome-relative representation
#' exome_sigs <- convert_signatures(cosmic_signatures_v2,
#'                                  opportunities_from = "human-genome",
#'                                  opportunities_to = "human-exome")
#' barplot(exome_sigs[1,])
#' 
#' # Normalise signatures to a sequence-independent representation
#' normalised_sigs <- convert_signatures(cosmic_signatures_v2,
#'                                       opportunities_from = "human-genome",
#'                                       opportunities_to = NULL)
#' barplot(normalised_sigs[1,])
#' 
#' # Convert back to genome-relative opportunities
#' genome_sigs <- convert_signatures(normalised_sigs,
#'                                   opportunities_from = NULL,
#'                                   opportunities_to = "human-genome")
#' barplot(genome_sigs[1,])
#' @export
convert_signatures <- function(signatures, opportunities_from = NULL, opportunities_to = NULL) {
    if (is.null(opportunities_from) & is.null(opportunities_to)) {
        stop("Either 'opportunities_from' or 'opportunities_to' must be provided.")
    }
    signatures <- to_matrix(signatures)
    opportunities_from <- build_opps_matrix(nrow(signatures), ncol(signatures), opportunities_from)
    opportunities_to <- build_opps_matrix(nrow(signatures), ncol(signatures), opportunities_to)
    
    conv_sigs <- signatures
    for (i in 1:nrow(conv_sigs)) {
        conv_sigs[i, ] <- conv_sigs[i, ] / opportunities_from[i, ] * opportunities_to[i, ]
        conv_sigs[i, ] <- conv_sigs[i, ] / sum(conv_sigs[i, ])
    }
    conv_sigs
}

#' Retrieve model parameters
#'
#' \code{retrieve_pars} obtains summary values for a set of model parameters (signatures, exposures,
#' activities or spectrum reconstructions) from a stanfit object.
#' @param mcmc_samples List with two elements named \code{`data`} and \code{`results`}, produced via
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}}, or
#' \code{\link{fit_extract_signatures}}.
#' @param par Name of the parameter set to extract. Admits character values \code{"signatures"},
#' \code{"exposures"}, \code{"activities"} or \code{"reconstructions"}.
#' @param hpd_prob Numeric value in the interval (0, 1), indicating the desired probability content
#' of HPD intervals (default is 0.95).
#' @return A list of three matrices, which contain the values corresponding to the means of the
#' model parameters and to the lower and upper ends of their HPD intervals, respectively.
#' @examples
#' \dontrun{
#' # Load example mutational catalogues
#' data("counts_21breast")
#'
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(counts_21breast, nsignatures = 2, model = "emu",
#'                               opportunities = "human-genome", iter = 800)
#'
#' # Retrieve signatures
#' signatures <- retrieve_pars(samples, "signatures")
#'
#' # Retrieve exposures and activities using custom HPD intervals
#' exposures <- retrieve_pars(samples, "exposures", hpd_prob = 0.9)
#' activities <- retrieve_pars(samples, "activities", hpd_prob = 0.975)
#'
#' # Retrieve reconstructed catalogues (reconstructions)
#' reconstructions <- retrieve_pars(samples, "reconstructions")
#'
#' # Plot signatures, reconstructions and mean exposures
#' plot_spectrum(signatures)
#' plot_spectrum(reconstructions)
#' barplot(t(exposures$mean))
#' }
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @export
retrieve_pars <- function(mcmc_samples, par, hpd_prob = 0.95) {
    if (!par %in% c("signatures", "exposures", "activities", "reconstructions")) {
        stop("'par' only admits the values \"signatures\", \"exposures\", \"activities\" or \"reconstructions\".")
    }
    if (par == "reconstructions") {
        p <- get_reconstructions(mcmc_samples)
        par_summ <- list(mean = apply(p$reconstructions, c(1, 3), sum),
                         lower = p$hpds[, 1, ],
                         upper = p$hpds[, 2, ])
        for (i in 1:length(par_summ)) {
            dimnames(par_summ[[i]]) <- dimnames(mcmc_samples$data$counts_real)
        }
    }
    else {
        p <- extract(mcmc_samples$result, pars = par)[[par]]
        if ("signatures" %in% names(mcmc_samples$data)) {
            signature_names <- rownames(mcmc_samples$data$signatures)
        }
        else {
            signature_names <- NULL
        }

        # Multi-sample case
        if (length(dim(p)) > 2) {
            # Assign dimension names
            if (par == "signatures") {
                if (dim(p)[3] %in% c(96, 192)) {
                    strand <- dim(p)[3] == 192  # strand bias indicator
                    names2 <- mut_types(strand)
                }
                else {
                    if (is.null(colnames(mcmc_samples$data$counts_real))) {
                        names2 <- paste("Mutation type", 1:ncol(mcmc_samples$data$counts_real))
                    }
                    else {
                        names2 <- colnames(mcmc_samples$data$counts_real)
                    }
                }
                if (is.null(signature_names)) {
                    LETTERLABELS <- letterwrap(dim(p)[2])
                    names1 <- paste("Signature", LETTERLABELS[1:dim(p)[2]])
                }
                else {
                    names1 <- signature_names
                }
            }
            else if (par %in% c("exposures", "activities")) {
                names1 <- rownames(mcmc_samples$data$counts_real)
                if (is.null(signature_names)) {
                    LETTERLABELS <- letterwrap(dim(p)[3])
                    names2 <- paste("Signature", LETTERLABELS[1:dim(p)[3]])
                }
                else {
                    names2 <- signature_names
                }
            }
            # for signatures: Signatures x Categories matrix
            # for exposures: Samples x Signatures matrix
            par_summ <- list(as.data.frame(matrix(NA, nrow = dim(p)[2], ncol = dim(p)[3], dimnames = list(names1, names2))),
                             as.data.frame(matrix(NA, nrow = dim(p)[2], ncol = dim(p)[3], dimnames = list(names1, names2))),
                             as.data.frame(matrix(NA, nrow = dim(p)[2], ncol = dim(p)[3], dimnames = list(names1, names2))))
            names(par_summ) <- c("mean", paste0(c("lower_", "upper_"), hpd_prob * 100))
            for (i in 1:nrow(par_summ[[1]])) {
                hpd <- HPDinterval(as.mcmc(p[,i,]), prob = hpd_prob)
                par_summ[[1]][i,] <- colMeans(p[,i,])
                par_summ[[2]][i,] <- hpd[,1]
                par_summ[[3]][i,] <- hpd[,2]
            }
        }

        # Single-sample case (only possible when fitting)
        else {
            names1 <- NULL
            if (is.null(signature_names)) {
                LETTERLABELS <- letterwrap(dim(p)[3])
                names2 <- paste("Signature", LETTERLABELS[1:dim(p)[3]])
            }
            else {
                names2 <- signature_names
            }
            par_summ <- list(as.data.frame(matrix(NA, nrow = 1, ncol = dim(p)[2], dimnames = list(names1, names2))),
                             as.data.frame(matrix(NA, nrow = 1, ncol = dim(p)[2], dimnames = list(names1, names2))),
                             as.data.frame(matrix(NA, nrow = 1, ncol = dim(p)[2], dimnames = list(names1, names2))))
            names(par_summ) <- c("mean", paste0(c("lower_", "upper_"), hpd_prob * 100))
            for (i in 1:nrow(par_summ[[1]])) {
                hpd <- HPDinterval(as.mcmc(p), prob = hpd_prob)
                par_summ[[1]][1,] <- colMeans(p)
                par_summ[[2]][1,] <- hpd[,1]
                par_summ[[3]][1,] <- hpd[,2]
            }
        }
    }
    par_summ
}

#' Generate posterior predictive check values from a model
#' @param mcmc_samples List with two elements named \code{`data`} and \code{`results`}, produced via
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}}, or
#' \code{\link{fit_extract_signatures}}.
#' @importFrom "stats" rmultinom rpois
#' @export
simulate_ppc <- function(mcmc_samples) {
    if (grepl("nmf", mcmc_samples$result@model_name)) {
        e <- extract(mcmc_samples$result, pars = "probs")
        ppc <- array(0, dim(e$probs))
        for (i in 1:dim(e$probs)[1]) {
            for (j in 1:dim(e$probs)[2]) {
                ppc[i, j, ] <- t(
                    rmultinom(1, sum(mcmc_samples$data$counts_real[j, ]), e$probs[i, j, ])
                )
            }
        }
    }
    else {
        e <- extract(mcmc_samples$result, pars = "expected_counts")
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
