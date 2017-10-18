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
#' \{build_catalogues} generates a set of mutational catalogues from a table containing 
#' the base change and trinucleotide context of each single-nucleotide variant in every sample.
#' @param variants Matrix with one row per single-nucleotide variant, and four or five columns:
#' \itemizer{
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
retrieve_pars <- function(mcmc_samples, feature, hpd_prob = 0.95, signature_names = NULL) {
    
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
    feat_summ
}

#' Plot mutational spectra
#' 
#' \code{plot_spectrum} generates plots of one or more spectra, which can be either mutational 
#' catalogues or mutational signatures. If the spectra contain values above 1, the values will be 
#' interpreted as mutation counts (as in a catalogue); otherwise, they will be interpreted as 
#' mutation probabilities (as in a signature). If multiple spectra are provided, one plot per 
#' spectrum is produced.
#' @param spectra Either a vector with one element for each of the 96 mutation types, or a matrix 
#' with 96 columns and one row per signature/catalogue, or a list of signatures as produced by 
#' \code{\link{retrieve_pars}}. In the latter case, HPD intervals will also be plotted. 
#' Row names will be adopted as the sample/signature names.
#' @param name Name to include in the plot title; useful when plotting a single spectrum.
#' @param pdf_path If provided, the plots will be output to a PDF file with this path. The PDF 
#' size and graphical parameters will be automatically set to appropriate values.
#' @param max_y Fixed maximum limit of the y-axis (if necessary).
#' @examples
#' # Load example mutational catalogues
#' data("counts_21breast")
#' 
#' # Plot catalogues
#' plot_spectrum(counts_21breast, pdf_path = "Catalogues.pdf")
#' 
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(counts_21breast, nsignatures = 2, method = "emu",
#'                               opportunities = "human-genome", iter = 800)
#' 
#' # Retrieve extracted signatures
#' sigs <- retrieve_pars(samples, "signatures")
#' 
#' # Plot signatures
#' plot_spectrum(sigs, pdf_path = "Signatures.pdf")
#' @export
plot_spectrum <- function(spectra, name = NULL, pdf_path = NULL, max_y = NULL) {
    # Fetch HPD interval values, if present
    if (is.list(spectra) & "mean" %in% names(spectra)) {
        spec <- spectra$mean
        lwr <- spectra$lower
        upr <- spectra$upper
    }
    else {
        spec <- spectra
        lwr <- NULL
        upr <- NULL
    }
    # Force spectrum to matrix (96 columns)
    spec <- to_matrix(spec)
    stopifnot(ncol(spec) %in% c(96, 192))
    
    NCAT <- ncol(spec)      # number of categories
    NSAMP <- nrow(spec)     # number of samples
    strand <- NCAT == 192   # strand bias indicator (logical)
    counts <- any(spec > 1) # count data indicator
    
    # Plot each spectrum
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    STRANDCOL <- c("deepskyblue3", "red3")
    BACKCOL <- c("#00BFFF33", "#00000033", "#EE2C2C33", "#C2C2C24D", "#A2CD5A4D", "#EEB4B44D")
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    BACKLIM <- c(0, 46.5, 93, 139.5, 186, 232.5, 279)
    
    if (!is.null(pdf_path)) {
        pdf(pdf_path, width = 24, height = 10.5)
        par(mar = c(9, 8, 6, 2.75))
    }
    
    # Standard spectrum (NCAT=96)
    if (!strand) {
        for (i in 1:NSAMP) {
            if (is.null(max_y)) {
                FACTOR <- 1.2
                samp_max_y <- max(0.05,
                                  ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
            }
            else {
                samp_max_y <- max_y
            }
            # Plot spectrum bars
            bars <- barplot(spec[i,], 
                            names.arg = mut_types(),
                            col = rep(COLORS, each = 16), border = "white",
                            yaxt = "n", ylim = c(0, samp_max_y), xlim = c(-1, 116),
                            cex = 1.3, cex.axis = 1.5, cex.lab = 1.7, 
                            las = 2, xaxs = "i", family = "mono")
            # Plot axis
            if (counts) {
                axis(side = 2, las = 2, cex.axis = 1.25)
                label <- "Mutations"
                n_text <- paste0(" (", sum(spec[i,]), " mutations)")
            }
            else {
                axis(side = 2, at = seq(0, samp_max_y, 0.05), las = 2, cex.axis = 1.25)
                label <- "Mutation probability"
                n_text <- ""
            }
            if (is.null(name)) {
                nme <- rownames(spec)[i]
            }
            else {
                nme <- name
            }
            if (NSAMP > 1) {
                num <- paste0(" #", i)
            }
            else {
                num <- ""
            }
            mtext(label, side = 2, cex = 1.7, line = 4.5)
            title(paste0("Mutational spectrum", num, "\n", nme, n_text), 
                  line = 1.5, cex.main = 2)
            # Plot HPD intervals
            if (!is.null(lwr)) {
                arrows(bars, spec[i,], bars, lwr[i,], angle = 90, 
                       length = 0.03, lwd = 1.5, col = "gray35")
                arrows(bars, spec[i,], bars, upr[i,], angle = 90, 
                       length = 0.03, lwd = 1.5, col = "gray35")
            }
            # Plot mutation type labels
            rect(xleft = XL, xright = XR, ybottom = 0.95 * samp_max_y, ytop = samp_max_y, 
                 col = COLORS, border = "white")
            text(x = (XL + XR) / 2, y = 0.91 * samp_max_y, labels = TYPES, cex = 2.25)
        }
    }
    
    # Strand-wise spectrum (NCAT=192)
    else {
        for (i in 1:NSAMP) {
            if (is.null(max_y)) {
                FACTOR <- 1.25
                samp_max_y <- max(0.05,
                                  ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
            }
            else {
                samp_max_y <- max_y
            }
            # Plot background panes and mutation type labels
            barplot(rbind(spec[i, 1:(NCAT/2)], spec[i, (NCAT/2+1):NCAT]), beside = TRUE, col = NA, border = NA,
                    space = c(0.1, 0.8), xaxs = "i", yaxt = "n", xaxt = "n", ylim = c(0, samp_max_y), xlim = c(-3, 280))
            for (j in 1:length(COLORS)) {
                rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0, 
                     ytop = samp_max_y, col = BACKCOL[j], border = "white")
                rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0.95 * samp_max_y, 
                     ytop = samp_max_y, col = COLORS[j], border = "white")
                text(x = (BACKLIM[j] + BACKLIM[j+1]) / 2, y = 0.91 * samp_max_y, 
                     labels = TYPES[j], cex = 2.25, col = "black")
            }
            # Plot spectrum bars
            bars <- barplot(rbind(spec[i, 1:(NCAT/2)], 
                                  spec[i, (NCAT/2+1):NCAT]), 
                            names.arg = mut_types(),
                            col = STRANDCOL, border = "white", beside = TRUE,
                            space = c(0.1, 0.8), yaxt = "n",
                            ylim = c(0, samp_max_y), xlim = c(-3, 270),
                            cex = 1.3, cex.axis = 1.5, cex.lab = 1.7,
                            las = 2, xaxs = "i", family = "mono", add = TRUE)
            # Plot legend and axis
            legend("topright", legend = c("Transcribed strand", "Untranscribed strand"), bty = "n",
                   fill = STRANDCOL, border = STRANDCOL, cex = 1.5, inset = c(0.018, 0.105))
            if (counts) {
                axis(side = 2, las = 2, cex.axis = 1.25)
                label <- "Mutations"
                n_text <- paste0(" (", sum(spec[i,]), " mutations)")
            }
            else {
                axis(side = 2, at = seq(0, samp_max_y, 0.05), las = 2, cex.axis = 1.25)
                label <- "Mutation probability"
                n_text <- ""
            }
            if (is.null(name)) {
                nme <- rownames(spec)[i]
            }
            else {
                nme <- name
            }
            if (NSAMP > 1) {
                num <- paste0(" #", i)
            }
            else {
                num <- ""
            }
            mtext(label, side = 2, cex = 1.7, line = 4.5)
            title(paste0("Mutational spectrum", num, "\n", nme, n_text), 
                  line = 1.5, cex.main = 2)
            # Plot HPD intervals
            if (!is.null(lwr)) {
                bars <- as.numeric(t(bars))
                arrows(bars, spec[i,], bars, lwr[i,], angle = 90, 
                       length = 0.03, lwd = 1.5, col = "gray35")
                arrows(bars, spec[i,], bars, upr[i,], angle = 90, 
                       length = 0.03, lwd = 1.5, col = "gray35")
            }
        }
    }
    
    if (!is.null(pdf_path)) {
        dev.off()
    }
}

#' Plot signature exposures
#' 
#' \code{plot_exposures} plots the distribution of signature exposures across the samples.
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param exposures Either a matrix of signature exposures, with one row per sample and one column 
#' per signature, or a list of exposures as produced by \code{\link{retrieve_pars}}. 
#' Only needed if \code{mcmc_samples} is not provided.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}. Only needed if \code{exposures} is not provided.
#' @param pdf_path If provided, the plots will be output to a PDF file with this path. The PDF 
#' size and graphical parameters will be automatically set to appropriate values.
#' @param signature_names Vector containing the names of the signatures. Only used when plotting
#' exposures obtained through signature fitting (not extraction).
#' @param thresh Minimum probability value that should be reached by the lower end of exposure HPD
#' intervals. The exposures for which the lower HPD bound is below this value will be colored in grey.
#' @param hpd_prob A value in the interval (0, 1), giving the target probability content of 
#' the HPD intervals.
#' @param horiz_labels If \code{TRUE}, sample name labels will be displayed horizontally in the
#' barplots.
#' @param legend_pos Character indicating the position of the legend in the exposures barplot. Admits values
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, 
#' \code{"topright"}, \code{"right"} and \code{"center"}.
#' @param sig_color_palette Character vector of color names or hexadecimal codes to use for each signature.
#' Must have at least as many elements as the number of signatures.
#' @export
plot_exposures <- function(counts, exposures = NULL, mcmc_samples = NULL, pdf_path = NULL,
                           signature_names = NULL, thresh = 0.01, hpd_prob = 0.95,
                           horiz_labels = FALSE, legend_pos = "topright", sig_color_palette = NULL) {
    if (is.null(exposures) & is.null(mcmc_samples)) {
        stop("Either 'exposures' (matrix or list) or 'mcmc_samples' (stanfit object) must be provided.")
    }
    if (!is.null(mcmc_samples)) {
        exposures <- retrieve_pars(mcmc_samples, "exposures", signature_names = signature_names,
                                   hpd_prob = hpd_prob)
        lwr <- exposures$lower
        upr <- exposures$upper
    }
    else if (is.list(exposures) & "mean" %in% names(exposures)) {
        lwr <- exposures$lower
        upr <- exposures$upper
    }
    else {
        lwr <- NULL
        upr <- NULL
    }
    if (!is.null(signature_names)) {
        for (i in 1:length(exposures)) {
            colnames(exposures[[i]]) <- signature_names
        }
    }
    exposures <- to_matrix(exposures)
    stopifnot(nrow(counts) == nrow(exposures))
    
    NSAMP <- nrow(counts)
    NSIG <- ncol(exposures)
    LETTERLABELS <- letterwrap(NSIG)
    
    if (is.null(rownames(counts))) {
        rownames(exposures) <- paste("Sample", 1:NSAMP)
    }
    else {
        rownames(exposures) <- rownames(counts)
    }
    if (is.null(colnames(exposures))) {
        colnames(exposures) <- paste("Signature", LETTERLABELS[1:NSIG])
    }
    
    if (is.null(sig_color_palette)) {
        sigcols <- default_sig_palette(NSIG)
    }
    else {
        sigcols <- sig_color_palette[1:NSIG]
    }
    
    if (!is.null(pdf_path)) {
        # PDF width increases with number of samples
        pdf(pdf_path, width = max(NSAMP * 0.13, 12), height = ifelse(NSAMP > 1, 12, 7))
        par(mar = c(6, 0, 4, 0), oma = c(1, 6, 1, 0))
    }
    
    par(lwd = 0.5)
    if (NSAMP > 1) {
        par(mfrow = c(3, 1))
    }
    
    # Obtain global (average) exposures
    exposures_global <- colMeans(exposures)
    if (!is.null(lwr)) {
        lwr_global <- colMeans(lwr)
        upr_global <- colMeans(upr)
        max_y <- max(upr_global)
    }
    else {
        max_y <- max(exposures_global)
    }
    
    # Plot global exposures
    colours <- rep("dodgerblue4", NSIG)
    if (!is.null(lwr)) {
        colours[lwr_global < thresh] <- "grey"
    }
    bars <- barplot(exposures_global, col = colours, border = NA, cex.names = 1.1, 
                    cex.main = 1.4, ylim = c(0, max_y), axes = F, las = ifelse(NSIG > 8, 2, 1),
                    main = "Global signature exposures across sample set")
    axis(side = 2, cex.axis = 1.1, las = 2, line = -2)
    mtext(text = "Mutation fraction", side = 2, line = 2.5)
    if (!is.null(lwr)) {
        arrows(bars, exposures_global, bars, lwr_global, 
               angle = 90, length = 0.03, lwd = 1.5, col = "gray50")
        arrows(bars, exposures_global, bars, upr_global, 
               angle = 90, length = 0.03, lwd = 1.5, col = "gray50")
    }
    
    # If >1 sample: plot exposures per sample
    if (NSAMP > 1) {
        par(mar = c(9, 0, 4, 0))
        
        # Obtain absolute exposures as mutation counts
        muts <- rowSums(counts)
        exposures_abs <- exposures * muts
        
        # Plot absolute exposures
        las = ifelse(horiz_labels, 1, 2)
        bars <- barplot(t(exposures_abs), col = sigcols, las = las, lwd = 0.25,
                        cex.names = 0.8, cex.main = 1.4, axes = FALSE,
                        main = "Signature exposures per sample (absolute)")
        axis(side = 2, cex.axis = 1.1, las = 2, line = -2)
        mtext(text = "Mutations", side = 2, line = 2.5)
        
        # Legend
        # expand legend box horizontally if there are lots of signatures
        LEGENDCOLS <- max(2, ceiling(NSIG / 10))
        legend(legend_pos, bty = "n", ncol = LEGENDCOLS, xpd = TRUE, inset = c(0.035, 0),
               fill = sigcols, border = sigcols, legend = colnames(exposures))
        
        # Plot relative exposures
        bars <- barplot(t(exposures), col = sigcols, las = las, lwd = 0.25,
                        cex.names = 0.8, cex.main = 1.4, axes = FALSE,
                        main = "Signature exposures per sample (relative)")
        axis(side = 2, cex.axis = 1.1, las = 2, line = -2)
        mtext(text = "Mutation fraction", side = 2, line = 2.5)
    }
    
    par(mfrow = c(1, 1), lwd = 1)
    if (!is.null(pdf_path)) {
        dev.off()
    }
}

#' Plot mutational spectrum reconstructions
#' 
#' \code{plot_reconstruction} plots the reconstructions of the original mutational catalogues using the 
#' signatures and/or exposures obtained through extraction or fitting. If provided with multiple catalogues, 
#' it generates one plot per catalogue. Fitting or extraction results can be provided either as a single 
#' stanfit object (generated via \code{\link{fit_signatures}} or \code{\link{extract_signatures}}), 
#' or as separate signatures and exposures matrices (or lists produced via \code{\link{retrieve_pars}}). 
#' Only the former option allows the incorporation of HPD intervals to the reconstructed catalogue.
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}. Only needed if \code{signatures} or \code{exposures}
#' are not provided.
#' @param signatures Either a matrix of mutational signatures, with one row per signature and one
#' column for each of the 96 mutation types, or a list of signatures generated via
#' \code{\link{retrieve_pars}}. Only needed if \code{mcmc_samples} is not provided.
#' @param exposures Either a matrix of signature exposures, with one row per sample and one column 
#' per signature, or a list of exposures as produced by \code{\link{retrieve_pars}}. 
#' Only needed if \code{mcmc_samples} is not provided.
#' @param opportunities If \code{signatures} and/or \code{exposures} were obtained by extracting or fitting
#' signatures using the "EMu" model (\code{method = "emu"}), these should be the same opportunities used 
#' for extraction/fitting. Admits values \code{"human-genome"} and \code{"human-exome"}.
#' @param pdf_path If provided, the plots will be output to a PDF file with this path. The PDF 
#' size and graphical parameters will be automatically set to appropriate values.
#' @param legend_pos Character indicating the position of the legend in the reconstructed spectrum. Admits values
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, 
#' \code{"topright"}, \code{"right"} and \code{"center"}.
#' @param sig_color_palette Character vector of color names or hexadecimal codes to use for each signature.
#' Must have at least as many elements as the number of signatures.
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
#' exposures <- retrieve_pars(samples, "exposures")
#' 
#' # Plot reconstructed catalogues using stanfit object
#' plot_reconstruction(counts_21breast, mcmc_samples = samples, opportunities = "human-genome",
#'                     pdf_path = "Reconstructions_1.pdf")
#'                     
#' # Plot reconstructed catalogues using retrieved signatures and exposures
#' plot_reconstruction(counts_21breast, signatures = signatures, exposures = exposures,
#'                     opportunities = "human-genome", pdf_path = "Reconstructions_2.pdf")
#' @importFrom "rstan" extract
#' @importFrom "coda" as.mcmc HPDinterval
#' @export
plot_reconstruction <- function(counts, mcmc_samples = NULL, signatures = NULL, 
                                exposures = NULL, opportunities = NULL, pdf_path = NULL,
                                legend_pos = "topright", sig_color_palette = NULL) {
    # Force counts to matrix
    counts <- to_matrix(counts)
    stopifnot(ncol(counts) %in% c(96, 192))
    
    NCAT <- ncol(counts)   # number of categories
    NSAMP <- nrow(counts)  # number of samples
    strand <- NCAT == 192  # strand bias indicator (logical)
    
    if (is.null(opportunities)) {
        opportunities <- matrix(1, nrow = NSAMP, ncol = NCAT)
    }
    else if (is.character(opportunities) & opportunities == "human-genome") {
        opportunities <- build_opps_matrix(NSAMP, "genome", strand)
    }
    else if (is.character(opportunities) & opportunities == "human-exome") {
        opportunities <- build_opps_matrix(NSAMP, "exome", strand)
    }
    else if (!is.matrix(opportunities)) {
        opportunities <- as.matrix(opportunities)
    }
    stopifnot(all(dim(opportunities) == dim(counts)))
    
    if (is.null(mcmc_samples) & (is.null(exposures) | is.null(signatures))) {
        stop("Either 'mcmc_samples' (stanfit object), or both 'signatures' and 'exposures' (matrices or lists), must be provided.")
    }
    
    cat("Building reconstructed catalogues...\n")
    
    # Case A: matrices given instead of MCMC samples
    if (is.null(mcmc_samples)) {
        # Force signatures and exposures to matrices
        signatures <- to_matrix(signatures)
        exposures <- to_matrix(exposures)
        NSIG <- nrow(signatures)
        
        stopifnot(ncol(signatures) == NCAT)
        stopifnot(nrow(exposures) == NSAMP)
        stopifnot(ncol(exposures) == NSIG)
        
        # Create reconstructed catalogues
        reconstructions <- array(NA, dim = c(NSAMP, NSIG, NCAT))
        for (i in 1:NSAMP) {
            rec <- exposures[i,] * signatures
            rec <- rec * opportunities[i,]
            rec <- rec / sum(rec)
            reconstructions[i,,] <- rec * sum(counts[i,])
        }
    }
    
    # Case B: MCMC samples given instead of matrices
    else {
        e <- extract(mcmc_samples)
        NREP <- dim(e$exposures)[1]
        stopifnot(NSAMP == dim(e$exposures)[2])

        # For fitting cases, signatures must be provided
        if (!("signatures" %in% names(e))) {
            if (is.null(signatures)) {
                stop("'mcmc_samples' contains signature fitting results: a signatures matrix must be provided via 'signatures'")
            }
            signatures <- to_matrix(signatures)
            
            # Reshape signatures as simulated MCMC samples
            e$signatures <- aperm(
                array(signatures, 
                      dim = c(nrow(signatures), NCAT, NREP)),
                c(3, 1, 2)
            )
        }
        NSIG <- dim(e$signatures)[2]
        
        # Obtain mean exposures (for legend)
        exposures <- t(apply(e$exposures, 2, colMeans))
        
        # Create reconstructed catalogues
        reconstructions <- array(NA, dim = c(NSAMP, NSIG, NCAT))
        hpds <- array(NA, dim = c(NSAMP, 2, NCAT))
        for (sample in 1:NSAMP) {
            
            # For EMu results
            if ("multiplier" %in% names(e)) {
                opp <- matrix(rep(as.matrix(opportunities[sample, ]), NSIG),
                              nrow = NSIG,
                              byrow = TRUE)
                
                arr <- aperm(
                    sapply(1:NREP, function(i) {
                        e$exposures[i, sample, ] * 
                            e$signatures[i, , ] * 
                            e$multiplier[i, sample] * 
                            opp
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
    }
    
    cat("Plotting reconstructions for each sample...\n")
    
    # Plotting
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    STRANDCOL <- c("deepskyblue3", "red3")
    BACKCOL <- c("#00BFFF33", "#00000033", "#EE2C2C33", "#C2C2C24D", "#A2CD5A4D", "#EEB4B44D")
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    BACKLIM <- c(0, 46.5, 93, 139.5, 186, 232.5, 279)
    
    if (is.null(sig_color_palette)) {
        sigcols <- default_sig_palette(NSIG)
    }
    else {
        sigcols <- sig_color_palette[1:NSIG]
    }
    
    if (is.null(rownames(counts))) {
        rownames(counts) <- paste("Sample", 1:NSAMP)
    }
    dimnames(reconstructions)[[3]] <- mut_types(strand)
    
    if (is.null(rownames(signatures))) {
        LETTERLABELS <- letterwrap(NSIG)
        sig_names <- paste("Signature", LETTERLABELS[1:NSIG])
    }
    else {
        sig_names <- rownames(signatures)
    }
    
    if (!is.null(pdf_path)) {
        stopifnot(is.character(pdf_path))
        pdf(pdf_path, width = 24, height = 18)
        par(mar = c(6, 8, 6, 2.75), oma = c(3, 0, 2, 0))
    }
    par(mfrow = c(2, 1))
    
    # Default spectrum (NCAT=96)
    if (!strand) {
        for (i in 1:NSAMP) {
            if (is.null(mcmc_samples)) {
                max_y <- max(c(counts[i, ], colSums(reconstructions[i, , ]))) * 1.1
            }
            else {
                max_y <- max(c(counts[i, ], hpds[i, , ])) * 1.1
            }
            
            # Plot original catalogue
            plot_spectrum(counts[i, ], name = rownames(counts)[i], max_y = max_y)
            
            # Plot catalogue reconstruction
            bars <- barplot(reconstructions[i, , ], 
                            col = sigcols, border = "white",
                            yaxt = "n", ylim = c(0, max_y), xlim = c(-1, 116),
                            cex = 1.3, cex.axis = 1.5, cex.lab = 1, las = 2, 
                            xaxs = "i", family = "mono")
            axis(side = 2, las = 2, cex.axis = 1.25)
            mtext("Mutations", side = 2, cex = 1.5, line = 4.5)
            title(paste0("Reconstructed mutational spectrum\nCosine similarity = ", 
                         round(cosine_sim(counts[i,], colSums(reconstructions[i, , ])), 3)),
                  line = 1, cex.main = 2)
            # HPD intervals
            if (!is.null(mcmc_samples)) {
                arrows(bars, colSums(reconstructions[i, , ]), 
                       bars, hpds[i, 1, ], 
                       angle = 90, length = 0.03, lwd = 1.5, col = "gray35")
                arrows(bars, colSums(reconstructions[i, , ]), 
                       bars, hpds[i, 2, ], 
                       angle = 90, length = 0.03, lwd = 1.5, col = "gray35")
            }
            # Mutation type labels
            rect(xleft = XL, xright = XR, ybottom = max_y * 0.95, ytop = max_y, 
                 col = COLORS, border = "white")
            text(x = (XL + XR) / 2, y = max_y * 0.9, labels = TYPES, cex = 2.25)
            # Legend
            legend(legend_pos, inset = c(0, 0.13), ncol = 2,
                   legend = paste0(sig_names, " (", round(exposures[i, ], 3), ")"), 
                   fill = sigcols, border = "white", cex = 1.5, bty = "n")
        }
    }
    
    # Strand-wise spectrum (NCAT=192)
    else {
        for (i in 1:NSAMP) {
            if (is.null(mcmc_samples)) {
                max_y <- max(c(counts[i, ], colSums(reconstructions[i, , ]))) * 1.1
            }
            else {
                max_y <- max(c(counts[i, ], hpds[i, , ])) * 1.1
            }
            
            # Plot original catalogue
            plot_spectrum(counts[i, ], name = rownames(counts)[i], max_y = max_y)
            
            # Plot catalogue reconstruction
            # Background panes and mutation type labels
            bars <- barplot(rbind(reconstructions[i, 1, 1:(NCAT/2)], 
                                  reconstructions[i, 1, (NCAT/2+1):NCAT]), 
                            names.arg = mut_types(), beside = TRUE, col = NA, border = NA, 
                            space = c(0.1, 0.8), xaxs = "i", yaxt = "n", ylim = c(0, max_y), 
                            xlim = c(-3, 280), family = "mono", cex = 1.3, cex.axis = 1.5, 
                            cex.lab = 1.7, las = 2)
            for (j in 1:length(COLORS)) {
                rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0, 
                     ytop = max_y, col = BACKCOL[j], border = "white")
                rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0.95 * max_y, 
                     ytop = max_y, col = COLORS[j], border = "white")
                text(x = (BACKLIM[j] + BACKLIM[j+1]) / 2, y = 0.91 * max_y, 
                     labels = TYPES[j], cex = 2.25, col = "black")
            }
            # Spectrum bars
            for (j in NSIG:1) {
                rec = colSums(to_matrix(reconstructions[i, 1:j, ]))
                barplot(rbind(rec[1:(NCAT/2)], rec[(NCAT/2+1):NCAT]), 
                        col = sigcols[j], border = "white", beside = TRUE,
                        space = c(0.1, 0.8), yaxt = "n", xaxt = "n",
                        ylim = c(0, max_y), xlim = c(-3, 270),
                        xaxs = "i", add = TRUE)
            }
            axis(side = 2, las = 2, cex.axis = 1.25)
            mtext("Mutations", side = 2, cex = 1.5, line = 4.5)
            title(paste0("Reconstructed mutational spectrum\nCosine similarity = ", 
                         round(cosine_sim(counts[i,], colSums(reconstructions[i, , ])), 3)),
                  line = 1, cex.main = 2)
            # HPD intervals
            if (!is.null(mcmc_samples)) {
                bars <- as.numeric(t(bars))
                arrows(bars, colSums(reconstructions[i, , ]), 
                       bars, hpds[i, 1, ], 
                       angle = 90, length = 0.03, lwd = 1.5, col = "gray35")
                arrows(bars, colSums(reconstructions[i, , ]), 
                       bars, hpds[i, 2, ], 
                       angle = 90, length = 0.03, lwd = 1.5, col = "gray35")
            }
            # Legend
            legend(legend_pos, inset = c(0.01, 0.105), ncol = 2,
                   legend = paste0(sig_names, " (", round(exposures[i, ], 3), ")"), 
                   fill = sigcols, border = sigcols, cex = 1.5, bty = "n")
        }
    }
    
    par(mfrow = c(1, 1))
    if (!is.null(pdf_path)) {
        dev.off()
    }
}

#' Plot all results from signature fitting or extraction
#' 
#' For a given set of signature fitting or extraction results, \code{plot_all} plots, in PDF format: 
#' \itemize{
#'  \item{All the original (input) mutational catalogues (via \code{\link{plot_spectrum}})}
#'  \item{Mutational signatures (via \code{\link{plot_spectrum}})}
#'  \item{Signature exposures (via \code{\link{plot_exposures}})}
#'  \item{All the reconstructed mutational spectra (via \code{\link{plot_reconstruction}})}
#' }
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param out_path Path to the directory where the output PDF files will be stored. M
#' Will be created if it does not exist.
#' @param prefix Optional prefix to be added to the output file names.
#' @param mcmc_samples Object of class stanfit, generated via either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}. Only needed if \code{signatures} or \code{exposures}
#' are not provided.
#' @param signatures Either a matrix of mutational signatures, with one row per signature and one
#' column for each of the 96 mutation types, or a list of signatures generated via
#' \code{\link{retrieve_pars}}. Only needed if \code{mcmc_samples} is not provided, or if it
#' contains signature fitting (instead of signature extraction) results.
#' @param exposures Either a matrix of signature exposures, with one row per sample and one column 
#' per signature, or a list of exposures as produced by \code{\link{retrieve_pars}}. 
#' Only needed if \code{mcmc_samples} is not provided.
#' @param opportunities If \code{signatures} and/or \code{exposures} were obtained by extracting or fitting
#' signatures using the "EMu" model (\code{method = "emu"}), these should be the same opportunities used 
#' for extraction/fitting. Admits values \code{"human-genome"} and \code{"human-exome"}.
#' @param thresh Minimum probability value that should be reached by the lower end of exposure HPD
#' intervals. The exposures for which the lower HPD bound is below this value will be colored in grey.
#' This value is passed to the \code{plot_exposures} function.
#' @param horiz_labels If \code{TRUE}, sample name labels will be displayed horizontally in the
#' barplots. This value is passed to the \code{plot_exposures} function.
#' @param hpd_prob A value in the interval (0, 1), giving the target probability content of 
#' the HPD intervals. This value is passed to the \code{plot_exposures} function.
#' @param exp_legend_pos Character indicating the position of the legend in the exposures barplot. Admits values
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, 
#' \code{"topright"}, \code{"right"} and \code{"center"}. This value is passed to the \code{plot_exposures} function.
#' @param rec_legend_pos Character indicating the position of the legend in the reconstructed spectrum. Admits values
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, 
#' \code{"topright"}, \code{"right"} and \code{"center"}. This value is passed to the \code{plot_reconstruction} function.
#' @param sig_color_palette Character vector of color names or hexadecimal codes to use for each signature.
#' Must have at least as many elements as the number of signatures. This value is passed to the 
#' \code{plot_exposures} and \code{plot_reconstruction} functions.
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
#' exposures <- retrieve_pars(samples, "exposures")
#' 
#' # Plot results using stanfit object
#' plot_all(counts_21breast, out_path = ".", prefix = "Test1", 
#'          mcmc_samples = samples, opportunities = "human-genome")
#'                     
#' # Plot results using retrieved signatures and exposures
#' plot_all(counts_21breast, out_path = ".", prefix = "Test2", 
#'          signatures = signatures, exposures = exposures, 
#'          opportunities = "human-genome")
#' @importFrom "rstan" extract
#' @export
plot_all <- function(counts, out_path, prefix = NULL, mcmc_samples = NULL, signatures = NULL, exposures = NULL,
                     opportunities = NULL, thresh = 0.01, horiz_labels = FALSE, hpd_prob = 0.95, 
                     signature_names = NULL, exp_legend_pos = "topright", rec_legend_pos = "topright", 
                     sig_color_palette = NULL) {
    
    if (is.null(mcmc_samples) & (is.null(exposures) | is.null(signatures))) {
        stop("Either 'mcmc_samples' (stanfit object), or both 'signatures' and 'exposures' (matrices or lists), must be provided.")
    }
    if (!is.null(mcmc_samples)) {
        if (is.null(signatures) & !("signatures" %in% mcmc_samples@model_pars)) { 
            stop("'mcmc_samples' contains signature fitting results: a signatures matrix must be provided via 'signatures'")
        }
    }
    
    # Create output directory if it does not exist
    dir.create(out_path, showWarnings=F)
    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    cat("Plotting original catalogues...\n")
    plot_spectrum(counts,
                  pdf_path = file.path(out_path, paste0(prefix, "Catalogues_", Sys.Date(), ".pdf")))
    
    # Case A: matrices provided instead of MCMC samples
    if (is.null(mcmc_samples)) {
        cat("Plotting mutational signatures...\n")
        plot_spectrum(signatures, 
                      pdf_path = file.path(out_path, paste0(prefix, "Signatures_", Sys.Date(), ".pdf")))
        
        cat("Plotting signature exposures...\n")
        plot_exposures(counts, exposures = exposures, signature_names = signature_names, 
                       thresh = thresh, sig_color_palette = sig_color_palette, 
                       horiz_labels = horiz_labels, legend_pos = exp_legend_pos,
                       pdf_path = file.path(out_path, paste0(prefix, "Exposures_", Sys.Date(), ".pdf")))
        
        plot_reconstruction(counts, signatures = signatures, exposures = exposures, 
                            opportunities = opportunities, legend_pos = rec_legend_pos, 
                            sig_color_palette = sig_color_palette,
                            pdf_path = file.path(out_path, paste0(prefix, "Reconstructions_", Sys.Date(), ".pdf")))
    }
    
    # Case B: MCMC samples provided instead of matrices
    else {
        cat("Plotting mutational signatures...\n")
        if ("signatures" %in% mcmc_samples@model_pars) {
            signatures <- retrieve_pars(mcmc_samples, feature = "signatures")
        }
        plot_spectrum(signatures, 
                      pdf_path = file.path(out_path, paste0(prefix, "Signatures_", Sys.Date(), ".pdf")))
        
        cat("Plotting signature exposures...\n")
        plot_exposures(counts, mcmc_samples = mcmc_samples, signature_names = signature_names, 
                       thresh = thresh, hpd_prob = hpd_prob, horiz_labels = horiz_labels, 
                       legend_pos = exp_legend_pos, sig_color_palette = sig_color_palette,
                       pdf_path = file.path(out_path, paste0(prefix, "Exposures_", Sys.Date(), ".pdf")))
        
        plot_reconstruction(counts, mcmc_samples = mcmc_samples, signatures = signatures, 
                            opportunities = opportunities, legend_pos = rec_legend_pos, 
                            sig_color_palette = sig_color_palette,
                            pdf_path = file.path(out_path, paste0(prefix, "Reconstructions_", Sys.Date(), ".pdf")))
    }
}

#' Plot goodness of fit
#' 
#' \code{plot_gof} plots the goodness of fit of a set of samples, each of which
#' has typically been sampled using an extraction model (EMu or NMF) with a 
#' different number of signatures.
#' @param sample_list List of objects of class stanfit. Elements which are not of class stanfit are ignored.
#' @param counts Matrix of observed mutation counts (integers), with one row per sample and 
#' column for each of the 96 mutation types.
#' @param stat Function for measuring goodness of fit. Admits values \code{"cosine"} 
#' (cosine similarity; default) or \code{"L2"} (L2 norm, a.k.a. Euclidean distance).
#' @importFrom "rstan" extract
#' @export
plot_gof <- function(sample_list, counts, stat = "cosine") {
    gof_function <- switch(stat,
                           "cosine" = cosine_sim,
                           "L2" = l2_norm)
    if (is.null(gof_function)) {
        stop("'stat' only admits values \"cosine\" and \"L2\".\nType ?plot_gof to see the documentation.")
    }
    
    nS <- NULL
    gof <- NULL
    for (samples in sample_list) {
        if (class(samples) != "stanfit") next
        
        if ("lambda" %in% samples@model_pars) {
            e <- extract(samples, pars = c("lambda", "signatures"))
            reconstructed <- apply(e$lambda, c(2, 3), mean)
            method <- "EMu"
        }
        else {
            e <- extract(samples, pars = c("probs", "signatures"))
            reconstructed <- apply(e$probs, c(2, 3), mean) * rowSums(counts)
            method = "NMF"
        }
        stopifnot(length(as.vector(reconstructed)) == length(as.vector(counts)))
        
        nS <- c(nS, dim(e$signatures)[2])
        gof <- c(gof, gof_function(as.vector(reconstructed),
                                        as.vector(counts)))
    }
    
    # Find the point of highest rate of change of gradient (i.e. highest positive 
    # 2nd derivative for 'L'-shaped curves, negative for 'r'-shaped curves)
    # Approximate 2nd derivative = x[n+1] + x[n-1] - 2x[n]
    deriv <- gof[3:(length(gof))] + gof[1:(length(gof)-2)] - 2 * gof[2:(length(gof)-1)]
    best <- ifelse(stat == "cosine", 
                   which.min(deriv) + 1,  # highest negative curvature
                   which.max(deriv) + 1)  # highest positive curvature
    
    plot(nS, gof, type = "o", lty = 3, pch = 16, col = "dodgerblue4",
         main = paste0("Goodness of fit (", stat, ")\nmodel: ", method),
         xlab = "Number of signatures", 
         ylab = paste0("Goodness of fit (", stat, ")"))
    points(nS[best], gof[best], pch = 16, col = "orangered", cex = 1.1)
    
    cat("Estimated best number of signatures:", nS[best], "\n")
    nS[best]
}

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
fit_signatures <- function(counts, signatures, exp_prior = NULL,
                           method = "nmf", opportunities = NULL, ...) {
    
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
    sampling(model, data = dat, ...)
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
#' @param nce NMF model specific: Use "normal-centred exposures" model (nce = TRUE). This parameterisation generally 
#' produces a higher ratio of effective sample size : chain length than the alternative, sampling exposures using a 
#' stick-breaking Dirichlet (used when nce = FALSE). Defaults to TRUE.
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
                               sig_prior = NULL, stanfunc = "sampling", nce = TRUE, ...) {
    
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
        
        if (nce) {
            model <- stanmodels$sigfit_ext_nmf_nce
        }
        else {
            model <- stanmodels$sigfit_ext_nmf
        }
        
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
                                   stanfunc = "sampling", ...) {
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
