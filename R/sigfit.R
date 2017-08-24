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

#' Returns default colour palette for signatures
default_sig_palette <- function(n) {
    c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
      "#A6761D", "#666666", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")[1:n]
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
        stop("`type` must be either \"genome\" or \"exome\"")
    }
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
    if (any(variants[,1] != sapply(strsplit(variants[,3], split=""), function(x) x[2]))) {
        stop("REF base (column 1) is not equal to middle base of the trinucleotide context (column 3).")
    }
    
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

#' Fetches COSMIC mutational signatures.
#' @param reorder Reorders the matrix by substitution type and trinucleotide.
#' @export
fetch_cosmic_data <- function(reorder = TRUE, remove_zeros = TRUE) {
    cosmic_sigs <- read.table("http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", 
                              header = TRUE, sep = "\t", check.names = FALSE)
    if (reorder) {
        cosmic_sigs <- cosmic_sigs[order(cosmic_sigs[["Substitution Type"]], cosmic_sigs[["Trinucleotide"]]),]
    }
    rownames(cosmic_sigs) <- cosmic_sigs[["Somatic Mutation Type"]]
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

#' Plot mutational spectra
#' Plots one or more spectra, which can be either mutational catalogues or mutational
#' signatures. If multiple spectra are provided, generates one plot per spectrum.
#' @param spectra Either a vector with 96 elements, a matrix with 96 columns and one row per signature/catalogue,
#' or a list of signature matrices, obtained via \code{retrieve_pars()}. In the latter case, HPD intervals will
#' also be plotted. Rownames are interpreted as the sample/signature names.
#' @param counts Do the values in \code{spectra} represent mutation counts instead of mutation probabilities?
#' @param name Name to include in the plot title; useful when plotting a single spectrum.
#' @examples
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "emu", 
#' opportunities = "human-genome")
#' 
#' # Retrieve and plot signatures
#' signatures <- retrieve_pars(samples, "signatures")
#' pdf("Signatures.pdf", width=24, height=11)
#' par(mar=c(9, 8, 6, 2.75))
#' plot_spectrum(signatures)
#' dev.off()
#' @export
plot_spectrum <- function(spectra, counts = FALSE, name = NULL, max_y = NULL, pdf_path = NULL) {
    NCAT <- 96  # number of categories
    # Fetch HPD interval values, if present
    if (is.list(spectra) & "mean" %in% names(spectra)) {
        spec <- spectra$mean
        lwr <- spectra$lower
        upr <- spectra$upper
    }
    else {
        spec <- spectra
        lwr <- NULL
    }
    # Force spectrum to matrix (96 columns)
    if (is.vector(spec)) {
        spec <- matrix(spec, nrow = 1)
    }
    if (!is.matrix(spec)) {
        spec <- as.matrix(spec)
    }
    stopifnot(ncol(spec) == NCAT)
    
    # Plot each spectrum
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    if (is.null(max_y)) {
        FACTOR <- ifelse(counts, 1.3, 1.1)
        max_y <- ifelse(is.null(upr), max(spec) * FACTOR, max(upr) * FACTOR)
    }
    if (!counts) {
        probs <- seq(0, 1, 0.05)
        max_y <- probs[which.max(probs > max_y)]
    }
    
    if (!is.null(pdf_path)) {
        pdf(pdf_path, width=24, height=11)
        par(mar=c(9, 8, 6, 2.75))
    }
    
    for (i in 1:nrow(spec)) {
        # Plot spectrum bars
        bars <- barplot(spec[i,], 
                        names.arg = mut_types(),
                        col = rep(COLORS, each = 16), border = "white",
                        yaxt = "n", ylim = c(0, max_y), xlim = c(-1, 116),
                        cex = 1.3, cex.axis = 1.5, cex.lab = 1.7, 
                        las = 2, xaxs = "i", family = "mono")
        if (counts) {
            axis(side = 2, las = 2, cex.axis = 1.25)
            label <- "Mutations"
            n_text <- paste0(" (N = ", sum(spec[i,]), " mutations)")
        }
        else {
            axis(side = 2, at = seq(0, max_y, 0.05), las = 2, cex.axis = 1.25)
            label <- "Mutation probability"
            n_text <- ""
        }
        if (is.null(name)) {
            nme <- rownames(spec)[i]
        }
        else {
            nme <- name
        }
        if (nrow(spec) > 1) {
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
                   length = 0.05, lwd = 1.5, col = "gray35")
            arrows(bars, spec[i,], bars, upr[i,], angle = 90, 
                   length = 0.05, lwd = 1.5, col = "gray35")
        }
        # Plot mutation type labels
        rect(xleft = XL, xright = XR, ybottom = max_y * 0.95, ytop = max_y, 
             col = COLORS, border = "white")
        text(x = (XL + XR) / 2, y = max_y * 0.9, labels = TYPES, cex = 2.25)
    }
    
    if (!is.null(pdf_path)) {
        dev.off()
    }
}

#' Plot signature exposures
#' \code{plot_exposures} produces barplots that show the distribution of
#' signatures exposures across catalogues.
#' @export
plot_exposures <- function(exposures, relative = TRUE, pdf_path = NULL, sig_color_palette = NULL) {
    if (is.list(exposures) & "mean" %in% names(exposures)) {
        exps <- exposures$mean
        lwr <- exposures$lower
        upr <- exposures$upper
    }
    else {
        exps <- as.matrix(exposures)
        lwr <- NULL
    }
    
    if (!is.null(pdf_path)) {
        pdf(pdf_path, 20, 12)
        par(mar=c(25, 4, 4, 2))
    }
    
    if (is.null(sig_color_palette)) {
        sigcols <- default_sig_palette(NSIG)
    }
    else {
        sigcols <- sig_color_palette[1:NSIG]
    }
    bars = barplot(t(exposures.refitted[,,"mean"])[3:1,], col=sigcols, las=2, space=0, main="Signature exposures")
    legend("topright", fill=rev(sig.cols), legend=rev(colnames(exposures)), xpd=T, bty="n", inset=c(-0.01,0))
    
    if (!is.null(pdf_path)) {
        dev.off()
    }
}

#' Plot mutational spectrum reconstructions
#' 
#' \code{plot_reconstruction} plots the reconstructions of the original mutational catalogues using the 
#' signatures and/or exposures obtained through extraction or fitting. If provided with multiple catalogues, it generates one
#' plot per catalogue. Data can be provided either as a single stanfit object resulting from extraction/fitting, 
#' or as separate objects for mutation counts, signatures and exposures. The latter option does not allow
#' plotting of error bars in the reconstructed catalogue.
#' @param counts Matrix of observed mutation counts, with 96 columns and
#' one row per sample.
#' @param mcmc_samples Object of class stanfit, produced by either \code{\link{fit_signatures}}
#' or \code{\link{extract_signatures}}. Needed if \code{counts}, \code{signatures} or \code{exposures}
#' are not provided.
#' @param signatures Either a matrix of mutational signatures, with 96 columns and one row per
#' signature, or a list of signatures as produced by \code{\link{retrieve_pars}}. Only needed if 
#' \code{mcmc_samples} is not provided.
#' @param exposures Either a matrix of signature exposures, with one row per sample and one column 
#' per signature, or a list of exposures as produced by \code{\link{retrieve_pars}}. Only needed if \code{mcmc_samples} is not provided.
#' @param opportunities If signatures or exposures were obtained using \code{method="emu"},
#' the opportunities used for extraction/fitting need to be provided here. Admits values \code{"human-genome"}
#' and \code{"human-exome"}. Only needed if \code{mcmc_samples} is not provided.
#' @param pdf_path If provided, the plots will be output to a PDF file with this path.
#' @param sig_color_palette Character vector of colour names or hexadecimal codes to use for each signature.
#' Must have at least as many elements as the number of signatures.
#' @examples
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "emu", 
#' opportunities = "human-genome")
#' 
#' # Retrieve signatures and exposures
#' signatures <- retrieve_pars(samples, "signatures")
#' exposures <- retrieve_pars(samples, "exposures")
#' 
#' # Plot reconstructed catalogues
#' pdf("Reconstruction.pdf", width = 24, height = 18)
#' par(mar = c(8, 8, 7, 2.75))
#' plot_reconstruction(mycounts, signatures, exposures, opportunities = "human-genome")
#' dev.off()
#' @importFrom "rstan" extract
#' @importFrom "coda" as.mcmc HPDinterval
#' @export
plot_reconstruction <- function(counts, mcmc_samples = NULL, signatures = NULL, exposures = NULL, 
                                opportunities = NULL, pdf_path = NULL, sig_color_palette = NULL) {
    
    NCAT <- 96             # number of categories (enforcing trinucleotides)
    NSAMP <- nrow(counts)  # number of samples
    
    # Force counts to matrix
    if (is.vector(counts)) {
        counts <- matrix(counts, nrow = 1)
    }
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }
    stopifnot(ncol(counts) == NCAT)
    
    if (is.null(opportunities)) {
        opportunities <- matrix(1, nrow = NSAMP, ncol = NCAT)
    }
    else if (opportunities == "human-genome") {
        opportunities <- matrix(rep(human_trinuc_freqs("genome"), NSAMP),
                                nrow = NSAMP, ncol = ncol(counts), byrow = TRUE)
    }
    else if (opportunities == "human-exome") {
        opportunities <- matrix(rep(human_trinuc_freqs("exome"), NSAMP),
                                nrow = NSAMP, ncol = ncol(counts), byrow = TRUE)
    }
    stopifnot(all(dim(opportunities) == dim(counts)))
    
    # Case A: matrices given instead of MCMC samples
    if (is.null(mcmc_samples)) {
        stopifnot(!(is.null(signatures) | is.null(exposures)))
        # Force signatures and exposures to matrices
        if (is.list(signatures) & "mean" %in% names(signatures)) {
            signatures <- signatures$mean
        }
        if (is.list(exposures) & "mean" %in% names(exposures)) {
            exposures <- exposures$mean
        }
        if (is.vector(signatures)) {
            signatures <- matrix(signatures, nrow = 1)
        }
        if (is.vector(exposures)) {
            exposures <- matrix(exposures, nrow = 1)
        }
        if (!is.matrix(signatures)) {
            signatures <- as.matrix(signatures)
        }
        if (!is.matrix(exposures)) {
            exposures <- as.matrix(exposures)
        }
        stopifnot(ncol(signatures) == NCAT)
        stopifnot(nrow(exposures) == NSAMP)
        stopifnot(ncol(exposures) == nrow(signatures))
        NSIG <- nrow(signatures)
        
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
        #total <- NULL        
        e <- extract(mcmc_samples)
        NREP <- dim(e$exposures)[1]
        stopifnot(NSAMP == dim(e$exposures)[2])

        # For fitting cases, signatures must be provided
        if (!("signatures" %in% names(e))) {
            if (is.null(signatures)) {
                stop("`mcmc_samples` contains signature fitting results: a signatures matrix must be provided via `signatures`")
            }
            if (is.list(signatures) & "mean" %in% names(signatures)) {
                signatures <- signatures$mean
            }
            if (is.vector(signatures)) {
                signatures <- matrix(signatures, nrow = 1)
            }
            if (!is.matrix(signatures)) {
                signatures <- as.matrix(signatures)
            }
            
            # Reshape signatures as simulated MCMC samples
            e$signatures <- aperm(
                array(signatures, 
                      dim = c(nrow(signatures), NCAT, NREP)),
                c(3, 1, 2)
            )
        }
        NSIG <- dim(e$signatures)[2]
        
        # Create reconstructed catalogues
        reconstructions <- array(NA, dim = c(NSAMP, NSIG, NCAT))
        hpds <- array(NA, dim = c(NSAMP, 2, NCAT))
        for (sample in 1:NSAMP) {
            
            # For EMu results
            if ("multiplier" %in% names(e)) {
                opp <- matrix(rep(as.matrix(opportunities[1, ]), NSIG),
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
            # if (is.null(total))
            #     total <- arr
            # else
            #     total <- total + arr

            reconstructions[sample, , ] <- apply(arr, c(2, 3), mean)
            hpds[sample, , ] <- t(HPDinterval(
                as.mcmc(apply(arr, c(1, 3), sum))
            ))
        }
        
        # Retrieve mean exposures
        exposures <- retrieve_pars(e, "exposures")$mean
    }
    
    # Plotting
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    
    if (is.null(sig_color_palette)) {
        sigcols <- default_sig_palette(NSIG)
    }
    else {
        sigcols <- sig_color_palette[1:NSIG]
    }
    
    if (!is.null(pdf_path)) {
        stopifnot(is.character(pdf_path))
        pdf(pdf_path, width = 24, height = 18)
        par(mar = c(8, 8, 7, 2.75))
    }
    par(mfrow = c(2, 1))
    
    if (is.null(rownames(counts))) {
        rownames(counts) <- paste("Sample", 1:nrow(counts))
    }
    dimnames(reconstructions)[[3]] <- mut_types()
    
    for (i in 1:NSAMP) {
        # Plot original catalogue
        plot_spectrum(counts[i,], counts = TRUE, name = rownames(counts)[i])
        
        # Plot catalogue reconstruction
        if (is.null(mcmc_samples)) {
            max_y <- max(colSums(reconstructions[i, , ])) * 1.3
        }
        else {
            max_y <- max(hpds[i, , ]) * 1.2
        }
        # Bars
        bars <- barplot(reconstructions[i, , ], 
                        col = sigcols, lwd = 0.3,
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
                   angle = 90, length = 0.05, lwd = 1.5, col = "gray35")
            arrows(bars, colSums(reconstructions[i, , ]), 
                   bars, hpds[i, 2, ], 
                   angle = 90, length = 0.05, lwd = 1.5, col = "gray35")
        }
        # Mutation type labels
        rect(xleft = XL, xright = XR, ybottom = max_y * 0.95, ytop = max_y, 
             col = COLORS, border = "white")
        text(x = (XL + XR) / 2, y = max_y * 0.9, labels = TYPES, cex = 2.25)
        # Legend
        if (is.null(rownames(signatures))) {
            sig_names <- paste("Signature", LETTERS[1:NSIG])
        }
        else {
            sig_names <- rownames(signatures)
        }
        legend("topright", inset = c(0, 0.13),
               legend = paste0(rev(sig_names), " (", round(exposures[i, ], 3), ")"), 
               fill = rev(sigcols),
               cex = 1.5, bty = "n")
    }
    
    par(mfrow = c(1, 1))
    if (!is.null(pdf_path)) {
        dev.off()
    }
}

#' Plot goodness of fit
#' 
#' \code{plot_gof} plots the goodness of fit of a set of samples, each of which
#' has typically been sampled using an extraction model (EMu or NMF) with a 
#' different number of signatures.
#' @param sample_list List of objects of class stanfit. Elements of the list
#' which are not of class stanfit are ignored.
#' @param counts Matrix of mutation counts, with 96 columns and one row per catalogue.
#' @param stat Function for measuring goodness of fit. Admits values "cosine" 
#' (cosine similarity; default) or "L2" (L2 norm, a.k.a. Euclidean distance).
#' @importFrom "rstan" extract
#' @export
plot_gof <- function(sample_list, counts, stat = "cosine") {
    gof_function <- switch(stat,
                           "cosine" = cosine_sim,
                           "L2" = l2_norm)
    if (is.null(gof_function)) {
        stop("Enter a valid option for `stat` -> \"cosine\", \"L2\"")
    }
    
    nS <- c()
    gof <- c()
    for (samples in sample_list) {
        if (class(samples) != "stanfit") next
        
        e <- extract(samples)
        if ("lambda" %in% names(e)) {
            reconstructed <- apply(e$lambda, c(2, 3), mean)
            method <- "EMu"
        }
        else {
            reconstructed <- apply(e$probs, c(2, 3), mean) * rowSums(counts)
            method = "NMF"
        }
        stopifnot(length(as.vector(reconstructed)) == length(as.vector(counts)))
        
        nS <- append(nS, dim(e$signatures)[2])
        gof <- append(gof, gof_function(as.vector(reconstructed),
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

#' Retrieve model parameters
#' 
#' Obtains summary values for a set of model parameters (signatures or exposures) from a stanfit object.
#' @param mcmc_samples An object of class stanfit, resulting from MCMC sampling through either \code{fit_signatures} or \code{extract_signatures}.
#' @param feature Name of the parameter set to extract; either \code{"signatures"} or \code{"exposures"}.
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
#' # Plot signatures and mean exposures
#' plot_spectrum(signatures)
#' barplot(exposures$mean)   ## or barplot(exposures[[1]])
#' @importFrom "rstan" extract
#' @importFrom "coda" HPDinterval
#' @importFrom "coda" as.mcmc
#' @export
retrieve_pars <- function(mcmc_samples, feature, prob = 0.95, signature_names = NULL) {
    feat <- extract(mcmc_samples, pars = feature)[[feature]]
    # Multi-sample case
    if (length(dim(feat)) > 2) {
        # Assign dimension names
        if (feature == "signatures") {
            names2 <- mut_types()
            if (is.null(signature_names)) {
                names1 <- paste("Signature", LETTERS[1:dim(feat)[2]])
            }
            else {
                if (dim(feat)[2] != length(signature_names)) {
                    stop("`signature_names` must have length equal to the number of signatures")
                }
                names1 <- signature_names
            }
        }
        else if (feature == "exposures") {
            names1 <- NULL
            if (is.null(signature_names)) {
                names2 <- paste("Signature", LETTERS[1:dim(feat)[3]])
            }
            else {
                if (dim(feat)[3] != length(signature_names))  {
                    stop("`signature_names` must have length equal to the number of signatures")
                }
                names2 <- signature_names
            }
        }
        # for signatures: Signatures x Categories matrix
        # for exposures: Samples x Signatures matrix
        feat_summ <- list(matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)),
                          matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)),
                          matrix(NA, nrow = dim(feat)[2], ncol = dim(feat)[3], dimnames = list(names1, names2)))
        names(feat_summ) <- c("mean", paste0(c("lower_", "upper_"), prob * 100))
        for (i in 1:nrow(feat_summ[[1]])) {
            hpd <- HPDinterval(as.mcmc(feat[,i,]), prob = prob)
            feat_summ[[1]][i,] <- colMeans(feat[,i,])
            feat_summ[[2]][i,] <- hpd[,1]
            feat_summ[[3]][i,] <- hpd[,2]
        }
    } 
    # Single-sample case (only possible when fitting)
    else {
        names1 <- NULL
        if (is.null(signature_names)) {
            names2 <- paste("Signature", LETTERS[1:dim(feat)[3]])
        }
        else {
            if (dim(feat)[3] != length(signature_names)) {
                stop("`signature_names` must have length equal to the number of signatures")
            }
            names2 <- signature_names
        }
        feat_summ <- list(matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)),
                          matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)),
                          matrix(NA, nrow = 1, ncol = dim(feat)[2], dimnames = list(names1, names2)))
        names(feat_summ) <- c("mean", paste0(c("lower_", "upper_"), prob * 100))
        for (i in 1:nrow(feat_summ[[1]])) {
            hpd <- HPDinterval(as.mcmc(feat), prob = prob)
            feat_summ[[1]][1,] <- colMeans(feat)
            feat_summ[[2]][1,] <- hpd[,1]
            feat_summ[[3]][1,] <- hpd[,2]
        }
    }
    feat_summ
}

#' Run MCMC to fit signatures and estimate exposures
#' @param counts Matrix of mutation counts per category (columns) per genome sample (rows).
#' @param signatures Matrix of mutational signatures (rows) to be fitted.
#' @param prior Vector with one element per signature, to be used as the Dirichlet prior in the sampling chain. 
#' Default prior is uniform (uninformative).
#' @param method Model to sample from; either \code{"nmf"} or \code{"emu"}.
#' @param opportunities Optional matrix of mutational opportunities for \code{"emu"} method; must have same 
#' dimension as \code{counts}. If equals to \code{"human-genome"} or \code{"human-exome"}, the reference h
#' uman genome/exome opportunities will be used for every sample.
#' @param ... Arguments to pass to \code{rstan::sampling}.
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
    stopifnot(length(prior) == nrow(signatures))
    
    # Force counts to matrix
    if (is.vector(counts)) {
        counts <- matrix(counts, nrow = 1)
    }
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NSIG, NCAT]
    stopifnot(ncol(counts) == ncol(signatures))
    stopifnot(length(prior) == nrow(signatures))
    
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
        if (!is.matrix(opportunities)) {
            opportunities <- as.matrix(opportunities)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        # NEED TO IMPLEMENT alpha
        dat <- list(
            C = ncol(counts),
            S = nrow(signatures),
            G = nrow(counts),
            signatures = signatures,
            counts = counts,
            opps = opportunities,
            alpha = prior
        )
        model <- stanmodels$sigfit_fit_emu
    }
    
    else {
        dat <- list(
            C = ncol(counts),
            S = nrow(signatures),
            G = nrow(counts),
            signatures = signatures,
            counts = counts,
            alpha = prior
        )
        model <- stanmodels$sigfit_fit_nmf
    }
    sampling(model, data = dat, ...)
}

#' Extract signatures from a set of mutation counts
#' 
#' @param counts Matrix of mutation counts for each sample (rows) in each category
#' (columns).
#' @param nsignatures Number (or range of numbers) of signatures to extract.
#' @param method Either \code{"emu"} (default) or \code{"nmf"}.
#' @param opportunities Optional matrix of mutational opportunities for \code{"emu"} method; 
#' must have same dimension as \code{counts}. If equals to \code{"human-genome"} or 
#' \code{"human-exome"}, the reference human genome/exome opportunities will be used for every sample.
#' @param exposures_prior Single numeric value to use for all the exposure priors; by default, 0.5 (i.e. Jeffreys prior).
#' @param stanfunc \code{"sampling"}|\code{"optimizing"}|\code{"vb"} Choice of rstan inference strategy. 
#' \code{"sampling"} is the full Bayesian MCMC approach, and is the default. \code{"optimizing"}
#' returns the Maximum a Posteriori (MAP) point estimates via numerical optimization.
#' \code{"vb"} uses Variational Bayes to approximate the full posterior.
#' @param ... Any other parameters to pass to rstan.
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @importFrom "rstan" extract
#' @export
extract_signatures <- function(counts, nsignatures, method = "emu", 
                               opportunities = NULL, exposures_prior = 0.5, 
                               stanfunc = "sampling", ...) {
    # Force counts to matrix
    if (is.vector(counts)) {
        counts <- matrix(counts, nrow = 1)
    }
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }
    
    # EMu model
    if (method == "emu") {
        # Build opportunities matrix
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
        if (!is.matrix(opportunities)) {
            opportunities <- as.matrix(opportunities)
        }
        stopifnot(all(dim(opportunities) == dim(counts)))
        
        model <- stanmodels$sigfit_ext_emu
        dat <- list(
            C = ncol(counts),
            G = nrow(counts),
            S = nsignatures[1],
            counts = counts,
            opps = opportunities
        )
    }
    
    # NMF model
    else if (method == "nmf") {
        if (!is.null(opportunities)) {
            warning("Using \"nmf\" model: `opportunities` will be ignored.")
        }
        
        model <- stanmodels$sigfit_ext_nmf
        dat <- list(
            C = ncol(counts),
            G = nrow(counts),
            S = as.integer(nsignatures[1]),
            counts = counts,
            exposures_prior_val = as.numeric(exposures_prior)
        )
    }
    else {
        stop("`method` must be either \"emu\" or \"nmf\".")
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
        
        names(out) <- paste0("nsignatures=", 1:length(out))
        
        # Plot goodness of fit and best number of signatures
        out$best <- plot_gof(out, counts)
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
#' @param stanfunc \code{"sampling"}|\code{"optimizing"}|\code{"vb"} Choice of rstan 
#' inference strategy. \code{"sampling"} is the full Bayesian MCMC approach, and is the 
#' default. \code{"optimizing"} returns the Maximum a Posteriori (MAP) point estimates 
#' via numerical optimization. \code{"vb"} uses Variational Bayes to approximate the 
#' full posterior.
#' @param ... Any other parameters to pass through to rstan.
#' @useDynLib sigfit, .registration = TRUE
#' @importFrom "rstan" sampling
#' @importFrom "rstan" optimizing
#' @importFrom "rstan" vb
#' @export
fit_extract_signatures <- function(counts, signatures, num_extra_sigs, 
                                   stanfunc = "sampling", ...) {
    # Force counts to matrix
    if (is.vector(counts)) {
        counts <- matrix(counts, nrow = 1)
    }
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }
    
    # Add pseudocounts to signatures
    signatures <- remove_zeros_(signatures)
    
    # Check dimensions are correct. Should be:
    # counts[NSAMPLES, NCAT], signatures[NSIG, NCAT]
    stopifnot(ncol(counts) == ncol(signatures))
    
    dat <- list(
        C = ncol(counts),
        S = nrow(signatures),
        G = nrow(counts),
        N = as.integer(num_extra_sigs),
        fixed_sigs = signatures,
        counts = counts
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
