#' Colour warnings for posterior predictive checks
#' yellow = observed value is in top or bottom 2.5% of ppc simulations
#' red = observed value is outside ppc simulations
make_colors <- function(c, c_ppc, sample = 1) {
    n_ppc_sims <- dim(c_ppc)[1]
    v <- vector("character", ncol(c))
    for (i in 1:ncol(c)) {
        prop <- sum(c[sample, i] < c_ppc[, sample, i]) / n_ppc_sims
        if (prop < 0.001 | prop > 0.999) {
            v[i] <- "firebrick"
        }
        else if (prop < 0.025 | prop > 0.975) {
            v[i] <- "yellow3"
        }
        else {
            v[i] <- "black"
        }
    }
    v
}

#' Plot result of posterior predictive check
#' black = observed value is not extreme compared to the simulated distribution
#' yellow = observed value is in top or bottom 2.5% of ppc simulations
#' red = observed value is outside ppc simulations
#' @param c Matrix of observed counts
#' @param c_ppc Simulated posterior predictive counts obtained from
#' MCMC samples using \code{extract(samples)$counts_ppc}
#' @export
plot_ppc <- function(c, c_ppc, sample = 1) {
    plot(c[sample, ], type = "n")
    n_ppc_sims <- dim(c_ppc)[1]
    for (i in sample(1:n_ppc_sims, 500)) {
        lines(c_ppc[i, sample, ], pch=20, cex=0.3,
              col = rgb(.09, .45, .80, 0.1))
    }
    colors <- make_colors(c, c_ppc, sample)
    lines(c[sample, ], type = 'p', lwd=2, col = colors, pch = 20)
    lines(c[sample, ], type = 'h', lend=1, lwd=1, col = colors)
    legend('topright', 
           legend=c("PPC distribution", 
                    "Observation (not extreme relative to PPC)", 
                    "Observation (in 5% tails of PPC)", 
                    "Observation (outside PPC)"), 
           col=c(rgb(.09, .45, .80, 1),
                 "black",
                 "yellow3",
                 "firebrick"), 
           lty=c(1, NA, NA, NA), 
           pch=c(NA, 20, 20, 20), 
           lwd=c(3, 2, 2, 2),
           cex=0.8)
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
#' @param counts Integer matrix of observed mutation counts, with one row per sample and 
#' column for each of the 96 mutation types.
#' @param stat Character, function for measuring goodness of fit. Admits values \code{"cosine"} 
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
        
        if (grepl("emu", samples@model_name)) {
            method <- "EMu"
            e <- extract(samples, pars = c("expected_counts", "signatures"))
            reconstructed <- apply(e$expected_counts, c(2, 3), mean)
        }
        
        else {
            method <- "NMF"
            e <- extract(samples, pars = c("probs", "signatures"))
            reconstructed <- apply(e$probs, c(2, 3), mean) * rowSums(counts)
        }
        
        stopifnot(dim(reconstructed) == dim(counts))
        
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
        if (!is.null(mcmc_samples) & grepl("emu", mcmc_samples@model_name)) {
            cat("Warning: Plotting EMu results, but no opportunities provided\n")
        }
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
        l <- get_reconstructions(counts, mcmc_samples, signatures)
        reconstructions <- l$reconstructions
        exposures <- l$exposures
        hpds <- l$hpds
        NSIG <- dim(l$exposures)[2]
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
#' @importFrom "graphics" arrows axis barplot legend lines mtext par plot points rect text title
#' @importFrom "grDevices" pdf dev.off rgb
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
