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
#' 
#' black = observed value is not extreme compared to the simulated distribution
#' yellow = observed value is in top or bottom 2.5% of ppc simulations
#' red = observed value is outside ppc simulations
#' @param c Integer matrix of observed counts.
#' @param c_ppc Simulated posterior predictive counts obtained from
#' MCMC samples using \code{extract(samples)$counts_ppc}.
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
    legend("topright",
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

#' Plot goodness-of-fit
#'
#' \code{plot_gof} plots the goodness-of-fit of a set of samples, each of which has typically been
#' sampled using an signature extraction model with a different number of signatures.
#' @param sample_list List containing the results of signature extraction using
#' \code{\link{extract_signatures}} with multiple numbers of signatures (see argument
#' \code{nsignatures} in \code{\link{extract_signatures}}).
#' @param stat Function for measuring goodness-of-fit. Admits character values \code{"cosine"}
#' (cosine similarity, default) or \code{"L2"} (L2 norm or Euclidean distance).
#' @importFrom "rstan" extract
#' @export
plot_gof <- function(sample_list, stat = "cosine") {
    
    if (sum(sapply(sample_list, is.list)) < 4) {
        warning("Goodness-of-fit analysis omitted when using less than 4 values of 'nsignatures'.")
    }
    
    else {
        gof_function <- switch(stat,
                               "cosine" = cosine_sim,
                               "L2" = l2_norm)
        if (is.null(gof_function)) {
            stop("'stat' only admits values \"cosine\" and \"L2\".\nType ?plot_gof to read the documentation.")
        }
    
        nS <- NULL
        gof <- NULL
        for (samples in sample_list) {
            if (!is.list(samples)) next
    
            counts <- samples$data$counts_real
    
            model <- samples$result@model_name
            e <- extract(samples$result, pars = c("expected_counts", "signatures"))
    
            if (samples$data$family == 1) {
                reconstructed <- apply(e$expected_counts, c(2, 3), mean) * rowSums(counts)
            }
            else {
                reconstructed <- apply(e$expected_counts, c(2, 3), mean)
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
             main = paste0("Goodness-of-fit (", stat, ")\nmodel: ", model),
             xlab = "Number of signatures",
             ylab = paste0("Goodness-of-fit (", stat, ")"))
        points(nS[best], gof[best], pch = 16, col = "orangered", cex = 1.1)
    
        cat("Estimated best number of signatures:", nS[best], "\n")
        nS[best]
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
#' @param mcmc_samples List with two elements named \code{`data`} and \code{`results`}, produced via
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}}, or
#' \code{\link{fit_extract_signatures}}. This is the preferred option for supplying data and
#' results, but can be replaced by the combination of arguments \code{counts}, \code{signatures},
#' \code{exposures} and \code{opportunities}.
#' @param out_path Character indicating the path to the directory where the output PDF files will
#' be stored. The directory will be created if it does not exist.
#' @param prefix Character indicating an optional prefix to be added to output file names.
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Only needed if \code{mcmc_samples} is not provided.
#' @param signatures Either a numeric matrix of mutational signatures, with one row per signature
#' and one column per mutation type, or a list of matrices generated via \code{\link{retrieve_pars}}.
#' Only needed if \code{mcmc_samples} is not provided.
#' @param exposures Either a numeric matrix of signature exposures, with one row per sample and one
#' column per signature, or a list of matrices generated via \code{\link{retrieve_pars}}. Only
#' needed if \code{mcmc_samples} is not provided.
#' @param opportunities Numeric matrix of mutational opportunities, with one row per signature and
#' one column per mutation type. It also admits character values \code{"human-genome"} or
#' \code{"human-exome"}, in which case the mutational opportunities of the reference human
#' genome/exome will be used. Only needed if \code{mcmc_samples} is not provided and opportunities
#' were used during signature extraction or fitting.
#' @param thresh Numeric indicating the minimum threshold for the lower HPD limits of signature
#' exposures (default is 0.01). Exposures with a lower HPD limit below this value will be shown in
#' grey. This value is passed to \code{\link{plot_exposures}}.
#' @param hpd_prob Numeric value in the interval (0, 1), indicating the desired probability content
#' of HPD intervals (default is 0.95). This value is passed to \code{\link{plot_exposures}}.
#' @param signature_names Character vector containing the name of each signature. Only needed if
#' \code{mcmc_samples} is not provided and the exposures were obtained via signature fitting
#' (rather than extraction). This value is passed to \code{\link{plot_exposures}}.
#' @param exp_margin_bottom Numeric indicating the bottom margin of the exposures barplots, in
#' inches (default is 10.5). This value is passed to \code{\link{plot_exposures}}.
#' @param exp_legend_pos Character indicating the position of the legend in the exposures barplot.
#' Admits values \code{"top"}, \code{"bottom"}, \code{"center"}, \code{"left"}, \code{"right"}, 
#' \code{"topleft"}, \code{"topright"}, \code{"bottomleft"} and \code{"bottomright"} (default is
#' \code{"topleft"}). This value is passed to \code{\link{plot_exposures}}.
#' @param exp_legend_cex Numeric indicating the relative size of the legend in the exposures
#' barplot (default is 2). This value is passed to \code{\link{plot_exposures}}.
#' @param exp_cex_names Numeric indicating the relative size of sample labels in the exposures
#' barplot (default is 1.9). This value is passed to \code{\link{plot_exposures}}.
#' @param rec_legend_pos Character indicating the position of the legend in the spectrum
#' reconstruction plots. Admits values \code{"top"}, \code{"bottom"}, \code{"center"}, \code{"left"},
#' \code{"right"}, \code{"topleft"}, \code{"topright"}, \code{"bottomleft"} and \code{"bottomright"}
#' (default is \code{"topleft"}). This value is passed to \code{\link{plot_reconstruction}}.
#' @param rec_legend_cex Numeric indicating the relative size of the legend in the reconstruction
#' plots (default is 2). This value is passed to \code{\link{plot_reconstruction}}.
#' @param sig_color_palette Character vector of custom color names or hexadecimal codes to use for
#' each signature in exposure and reconstruction plots. Must have at least as many elements as the
#' number of signatures. This value is passed to \code{\link{plot_exposures}} and
#' \code{\link{plot_reconstruction}}.
#' @param boxes Logical indicating whether boxes should be drawn around spectrum, signature and
#' reconstruction plots (default is \code{TRUE}). This value is passed to
#' \code{\link{plot_spectrum}} and \code{\link{plot_reconstruction}}.
#' @examples
#' \dontrun{
#' # Load example mutational catalogues
#' data("counts_21breast")
#'
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(counts_21breast, nsignatures = 2, model = "emu",
#'                               opportunities = "human-genome", iter = 800)
#'
#' # Retrieve signatures and exposures
#' signatures <- retrieve_pars(samples, "signatures")
#' exposures <- retrieve_pars(samples, "exposures")
#'
#' # Plot results using MCMC samples
#' plot_all(mcmc_samples = samples, out_path = ".", prefix = "Test1")
#'
#' # Plot results using retrieved signatures and exposures
#' plot_all(counts = counts_21breast, signatures = signatures,
#'          exposures = exposures, opportunities = "human-genome",
#'          out_path = ".", prefix = "Test2")
#' }
#' @importFrom "rstan" extract
#' @export
plot_all <- function(mcmc_samples = NULL, out_path, prefix = NULL, counts = NULL, signatures = NULL,
                     exposures = NULL, opportunities = NULL, thresh = 0.01, hpd_prob = 0.95, 
                     signature_names = NULL, exp_margin_bottom = 10.5, exp_legend_pos = "topleft",
                     exp_legend_cex = 2, exp_cex_names = 1.9, rec_legend_pos = "topleft",
                     rec_legend_cex = 2, sig_color_palette = NULL, boxes = TRUE) {

    if (is.null(mcmc_samples) &
        (is.null(counts) | is.null(signatures) | is.null(exposures))) {
        stop("Either 'mcmc_samples', or all of 'counts', 'signatures' and 'exposures', must be provided.")
    }

    # Create output directory if it does not exist
    dir.create(out_path, showWarnings=F)
    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "_")
    }

    # Case A: matrices provided instead of MCMC samples
    if (is.null(mcmc_samples)) {
        cat("Plotting original catalogues...\n")
        plot_spectrum(counts, boxes = boxes,
                      pdf_path = file.path(out_path, paste0(prefix, "Catalogues_", Sys.Date(), ".pdf")))

        cat("Plotting mutational signatures...\n")
        plot_spectrum(signatures, boxes = boxes,
                      pdf_path = file.path(out_path, paste0(prefix, "Signatures_", Sys.Date(), ".pdf")))

        cat("Plotting signature exposures...\n")
        plot_exposures(counts = counts, exposures = exposures,
                       signature_names = signature_names, thresh = thresh,
                       sig_color_palette = sig_color_palette, cex_names = exp_cex_names,
                       margin_bottom = exp_margin_bottom, legend_pos = exp_legend_pos,
                       legend_cex = exp_legend_cex,
                       pdf_path = file.path(out_path, paste0(prefix, "Exposures_", Sys.Date(), ".pdf")))

        plot_reconstruction(counts = counts, signatures = signatures, exposures = exposures,
                            opportunities = opportunities, legend_pos = rec_legend_pos,
                            legend_cex = rec_legend_cex, sig_color_palette = sig_color_palette,
                            pdf_path = file.path(out_path, paste0(prefix, "Reconstructions_", Sys.Date(), ".pdf")))
    }

    # Case B: MCMC samples provided instead of matrices
    else {
        cat("Plotting original catalogues...\n")
        plot_spectrum(mcmc_samples$data$counts_real, boxes = boxes,
                      pdf_path = file.path(out_path, paste0(prefix, "Catalogues_", Sys.Date(), ".pdf")))

        cat("Plotting mutational signatures...\n")
        if ("signatures" %in% mcmc_samples$result@model_pars) {
            signatures <- retrieve_pars(mcmc_samples, "signatures")
        }
        else {
            signatures <- mcmc_samples$data$signatures
        }
        plot_spectrum(signatures, boxes = boxes,
                      pdf_path = file.path(out_path, paste0(prefix, "Signatures_", Sys.Date(), ".pdf")))

        cat("Plotting signature exposures...\n")
        plot_exposures(mcmc_samples = mcmc_samples,
                       thresh = thresh, hpd_prob = hpd_prob,
                       margin_bottom = exp_margin_bottom, legend_pos = exp_legend_pos,
                       sig_color_palette = sig_color_palette, cex_names = exp_cex_names,
                       legend_cex = exp_legend_cex,
                       pdf_path = file.path(out_path, paste0(prefix, "Exposures_", Sys.Date(), ".pdf")))

        plot_reconstruction(mcmc_samples = mcmc_samples,
                            legend_pos = rec_legend_pos, legend_cex = rec_legend_cex,
                            sig_color_palette = sig_color_palette, boxes=boxes,
                            pdf_path = file.path(out_path, paste0(prefix, "Reconstructions_", Sys.Date(), ".pdf")))
    }
}

#' Plot mutational spectrum reconstructions
#'
#' \code{plot_reconstruction} plots reconstructions of the original mutational catalogues obtained
#' using the inferred signatures and/or exposures. If provided with multiple catalogues, it produces
#' one plot per catalogue. Fitting or extraction results can be provided either as a single stanfit
#' object (generated via \code{\link{fit_signatures}} or \code{\link{extract_signatures}}), or as
#' separate signatures and exposures matrices (or lists produced via \code{\link{retrieve_pars}}).
#' Only the former option allows the incorporation of HPD intervals to the reconstructed catalogue.
#' @param mcmc_samples List with two elements named \code{`data`} and \code{`results`}, produced via
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}}, or
#' \code{\link{fit_extract_signatures}}. This is the preferred option for supplying data and
#' results, but can be replaced by the combination of arguments \code{counts}, \code{signatures},
#' \code{exposures} and \code{opportunities}.
#' @param pdf_path Character indicating the path to an optional output PDF file for the plots. The
#' PDF dimensions and graphical parameters are automatically set to appropriate values, but custom
#' dimensions can be specified via the arguments \code{pdf_width} and \code{pdf_height}.
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Only needed if \code{mcmc_samples} is not provided.
#' @param signatures Either a numeric matrix of mutational signatures, with one row per signature
#' and one column per mutation type, or a list of matrices generated via \code{\link{retrieve_pars}}.
#' Only needed if \code{mcmc_samples} is not provided.
#' @param exposures Either a numeric matrix of signature exposures, with one row per sample and one
#' column per signature, or a list of matrices generated via \code{\link{retrieve_pars}}. Only
#' needed if \code{mcmc_samples} is not provided.
#' @param opportunities Numeric matrix of mutational opportunities, with one row per signature and
#' one column per mutation type. It also admits character values \code{"human-genome"} or
#' \code{"human-exome"}, in which case the mutational opportunities of the reference human
#' genome/exome will be used. Only needed if \code{mcmc_samples} is not provided and opportunities
#' were used during signature extraction or fitting.
#' @param pdf_width Numeric indicating the width of the output PDF, in inches (default is 24).
#' Only used if \code{pdf_path} is provided.
#' @param pdf_height Numeric indicating the height of the output PDF, in inches (default is 15).
#' Only used if \code{pdf_path} is provided.
#' @param legend_pos Character indicating the position of the legend in the plots. Admits values
#' \code{"top"}, \code{"bottom"}, \code{"center"}, \code{"left"}, \code{"right"}, \code{"topleft"},
#' \code{"topright"}, \code{"bottomleft"} and \code{"bottomright"} (default is \code{"topleft"}).
#' @param legend_cex Numeric indicating the relative size of the legend (default is 2).
#' @param sig_color_palette Character vector of custom color names or hexadecimal codes to use for
#' each signature in exposure and reconstruction plots. Must have at least as many elements as the
#' number of signatures.
#' @param boxes Logical indicating whether boxes should be drawn around the plots (default is
#' \code{TRUE}).
#' @examples
#' \dontrun{
#' # Load example mutational catalogues
#' data("counts_21breast")
#'
#' # Extract signatures using the EMu (Poisson) model
#' samples <- extract_signatures(counts_21breast, nsignatures = 2, model = "emu",
#'                               opportunities = "human-genome", iter = 800)
#'
#' # Retrieve signatures and exposures
#' signatures <- retrieve_pars(samples, "signatures")
#' exposures <- retrieve_pars(samples, "exposures")
#'
#' # Plot reconstructed catalogues using stanfit object
#' plot_reconstruction(mcmc_samples = samples, pdf_path = "Reconstructions_1.pdf")
#'
#' # Plot reconstructed catalogues using retrieved signatures and exposures
#' plot_reconstruction(counts = counts_21breast, signatures = signatures, exposures = exposures,
#'                     opportunities = "human-genome", pdf_path = "Reconstructions_2.pdf")
#' }
#' @importFrom "rstan" extract
#' @importFrom "coda" as.mcmc HPDinterval
#' @importFrom "grDevices" cairo_pdf
#' @export
plot_reconstruction <- function(mcmc_samples = NULL, pdf_path = NULL, counts = NULL,
                                signatures = NULL, exposures = NULL, opportunities = NULL, 
                                pdf_width = 24, pdf_height = 15, legend_pos = "topleft", 
                                legend_cex = 2, sig_color_palette = NULL, boxes = TRUE) {

    if (!is.null(mcmc_samples)) {
        counts <- mcmc_samples$data$counts_real
        opportunities <- mcmc_samples$data$opportunities
    }
    else {
        if (is.null(counts) | is.null(exposures) | is.null(signatures)) {
            stop("Either 'mcmc_samples', or all three of 'counts', 'signatures' and 'exposures', must be provided.")
        }
    }

    # Force counts to (integer) matrix
    counts <- to_matrix(counts, int = TRUE)

    NCAT <- ncol(counts)   # number of categories
    NSAMP <- nrow(counts)  # number of samples
    strand <- NCAT == 192  # strand bias indicator (logical)

    if (is.null(opportunities) | is.character(opportunities)) {
        opportunities <- build_opps_matrix(NSAMP, NCAT, opportunities)
    }
    else if (!is.matrix(opportunities)) {
        opportunities <- as.matrix(opportunities)
    }
    stopifnot(all(dim(opportunities) == dim(counts)))

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
        l <- get_reconstructions(mcmc_samples)
        reconstructions <- l$reconstructions
        exposures <- l$exposures
        hpds <- l$hpds
        NSIG <- dim(l$exposures)[2]
    }

    cat("Plotting reconstructed catalogues...\n")

    # Plotting
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    STRANDCOL <- c("deepskyblue3", "red3")
    BACKCOL <- c("#00BFFF33", "#00000033", "#EE2C2C33", "#C2C2C24D", "#A2CD5A4D", "#EEB4B44D")
    LINECOL <- "gray60"
    FACTOR <- 1.095
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

    if ("signatures" %in% names(mcmc_samples$data)) {
        sig_names <- rownames(mcmc_samples$data$signatures)
    }
    else if (!is.null(rownames(signatures))) {
        sig_names <- rownames(signatures)
    }
    else {
        LETTERLABELS <- letterwrap(NSIG)
        sig_names <- paste("Signature", LETTERLABELS[1:NSIG])
    }

    if (!is.null(pdf_path)) {
        stopifnot(is.character(pdf_path))
        cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, onefile = TRUE)
        par(oma = c(1, 0, 1, 0))
        if (ncol(counts) %in% c(96, 192)) {
            par(mar = c(4.5, 7, 6.5, 2))
        }
        else {
            par(mar = c(9, 7, 6.5, 2))
        }
    }
    par(mfrow = c(2, 1))

    # Generic spectrum (NCAT!={96,192})
    if (!(ncol(counts) %in% c(96, 192))) {
        if (is.null(colnames(counts))) {
            types <- paste("Mut. type", 1:ncol(counts))
        }
        else {
            types <- colnames(counts)
        }

        for (i in 1:NSAMP) {
            if (is.null(mcmc_samples)) {
                max_y <- max(c(counts[i, ], colSums(reconstructions[i, , ]))) * FACTOR
            }
            else {
                max_y <- max(c(counts[i, ], hpds[i, , ])) * FACTOR
            }

            # Plot original catalogue
            plot_spectrum(counts[i, ], name = rownames(counts)[i], max_y = max_y, boxes = boxes)

            # Plot catalogue reconstruction
            bars <- barplot(reconstructions[i, , ],
                            names.arg = types, col = sigcols, border = "white",
                            yaxt = "n", ylim = c(0, max_y), cex.names = 1, las = 2)
            axis(side = 2, cex.axis = 1.9, lwd = 2)
            mtext("Mutations", side = 2, cex = 2, line = 3.5)
            title(paste0("Reconstructed spectrum (cosine similarity = ",
                         round(cosine_sim(counts[i,], colSums(reconstructions[i, , ])), 3), ")"),
                  line = 4, cex.main = 2.5)
            # HPD intervals
            if (!is.null(mcmc_samples)) {
                arrows(bars, hpds[i, 1, ],
                       bars, hpds[i, 2, ],
                       length = 0, lwd = 3, col = LINECOL)
            }
            # Legend
            legend(legend_pos, inset = c(0, 0.13), ncol = 2,
                   legend = paste0(sig_names, " (", round(exposures[i, ], 3), ")"),
                   fill = sigcols, border = "white", cex = legend_cex, bty = "n")
        }
    }
    else {
        dimnames(reconstructions)[[3]] <- mut_types(strand)

        # Default spectrum (NCAT=96)
        if (!strand) {
            for (i in 1:NSAMP) {
                if (is.null(mcmc_samples)) {
                    max_y <- max(c(counts[i, ], colSums(reconstructions[i, , ]))) * FACTOR
                }
                else {
                    max_y <- max(c(counts[i, ], hpds[i, , ])) * FACTOR
                }
                if (boxes) {
                    xlim <- c(-0.105, 115.5)
                }
                else {
                    xlim <- c(-1, 116)
                }

                # Plot original catalogue
                plot_spectrum(counts[i, ], name = rownames(counts)[i], max_y = max_y, boxes = boxes)

                # Plot catalogue reconstruction
                bars <- barplot(reconstructions[i, , ],
                                names.arg = substr(mut_types(), 1, 3), mgp = c(3, 0.8, 0),
                                col = sigcols, border = "white",
                                yaxt = "n", ylim = c(0, max_y), xlim = xlim,
                                las = 2, cex.names = 1.6, xaxs = "i", family = "mono")
                for (j in 1:length(COLORS)) {
                    idx <- ((j-1) * 16 + 1):(j * 16)
                    axis(side = 1, at = bars[idx], tick = FALSE, cex.axis = 1.6,
                         mgp = c(3, 0.8, 0), las = 2, family = "mono", font = 2,
                         col.axis = COLORS[j], labels = paste0(" ", substr(mut_types()[idx], 2, 2), " "))
                }
                axis(side = 2, cex.axis = 1.9, lwd = 2)
                mtext("Mutations", side = 2, cex = 2.4, line = 3.5)
                title(paste0("Reconstructed spectrum (cosine similarity = ",
                             round(cosine_sim(counts[i,], colSums(reconstructions[i, , ])), 3), ")"),
                      line = 4, cex.main = 2.5)
                # HPD intervals
                if (!is.null(mcmc_samples)) {
                    arrows(bars, hpds[i, 1, ],
                           bars, hpds[i, 2, ],
                           length = 0, lwd = 3, col = LINECOL)
                }
                # Mutation type labels
                text(x = (XL + XR) / 2, y = max_y * 1.055,
                     labels = TYPES, cex = 2.4, xpd = TRUE)
                rect(xleft = XL, xright = XR, ybottom = max_y * 0.945,
                     ytop = max_y, col = COLORS, border = "white")
                # Legend
                legend(legend_pos, inset = c(0, 0.05), ncol = 2,
                       legend = paste0(sig_names, " (", round(exposures[i, ], 3), ")"),
                       fill = sigcols, border = sigcols, cex = 1.9, bty = "n")
                # Box
                if (boxes) {
                    box(lwd = 2)
                }
            }
        }

        # Strand-wise spectrum (NCAT=192)
        else {
            for (i in 1:NSAMP) {
                if (is.null(mcmc_samples)) {
                    max_y <- max(c(counts[i, ], colSums(reconstructions[i, , ]))) * FACTOR
                }
                else {
                    max_y <- max(c(counts[i, ], hpds[i, , ])) * FACTOR
                }
                if (boxes) {
                    xlim <- c(0, 279.2)
                }
                else {
                    xlim <- c(-3, 280)
                }

                # Plot original catalogue
                plot_spectrum(counts[i, ], name = rownames(counts)[i], max_y = max_y, boxes = boxes)
                
                # Plot catalogue reconstruction
                # Background panes and mutation type labels
                barplot(rbind(reconstructions[i, 1, 1:(NCAT/2)],
                              reconstructions[i, 1, (NCAT/2+1):NCAT]),
                        beside = TRUE, col = NA, border = NA,
                        space = c(0.1, 0.8), xaxs = "i", yaxt = "n", xaxt = "n", 
                        ylim = c(0, max_y), xlim = xlim)
                for (j in 1:length(COLORS)) {
                    rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0,
                         ytop = max_y, col = BACKCOL[j], border = "white")
                    rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0.945 * max_y,
                         ytop = max_y, col = COLORS[j], border = "white")
                    text(x = (BACKLIM[j] + BACKLIM[j+1]) / 2, y = 1.055 * max_y,
                         labels = TYPES[j], cex = 2.4, xpd = TRUE)
                }
                # Spectrum bars
                for (j in NSIG:1) {
                    rec = colSums(to_matrix(reconstructions[i, 1:j, ]))
                    bars <- barplot(rbind(rec[1:(NCAT/2)],
                                          rec[(NCAT/2+1):NCAT]),
                                    names.arg = substr(mut_types(), 1, 3),
                                    col = sigcols[j], border = "white", las = 2,
                                    beside = TRUE, space = c(0.1, 0.8), mgp = c(3, 0.8, 0), 
                                    cex.names = 1.6, yaxt = "n", xaxs = "i",
                                    family = "mono", add = TRUE)
                }
                for (j in 1:length(COLORS)) {
                    idx <- ((j-1) * 16 + 1):(j * 16)
                    axis(side = 1, tick = FALSE, at = colMeans(bars[, idx]),
                         cex.axis = 1.6, mgp = c(3, 0.8, 0), las = 2,
                         family = "mono", font = 2, col.axis = COLORS[j],
                         labels = paste0(" ", substr(mut_types()[idx], 2, 2), " "))
                }
                axis(side = 2, cex.axis = 1.9, lwd = 2)
                mtext("Mutations", side = 2, cex = 2.4, line = 3.5)
                title(paste0("Reconstructed spectrum (cosine similarity = ",
                             round(cosine_sim(counts[i,], colSums(reconstructions[i, , ])), 3), ")"),
                      line = 4, cex.main = 2.5)
                # HPD intervals
                if (!is.null(mcmc_samples)) {
                    bars <- as.numeric(t(bars))
                    arrows(bars, hpds[i, 1, ],
                           bars, hpds[i, 2, ],
                           length = 0, lwd = 2.5, col = LINECOL)
                }
                # Legend
                legend(legend_pos, inset = c(0, 0.05), ncol = 2,
                       legend = paste0(sig_names, " (", round(exposures[i, ], 3), ")"),
                       fill = sigcols, border = sigcols, cex = 1.9, bty = "n")
                # legend(legend_pos, inset = c(0.01, 0.105), ncol = 2,
                #        legend = paste0(sig_names, " (", round(exposures[i, ], 3), ")"),
                #        fill = sigcols, border = sigcols, cex = 1.8, bty = "n")
                # Box
                if (boxes) {
                    box(lwd = 2)
                }
            }
        }
    }
    par(mfrow = c(1, 1))
    if (!is.null(pdf_path)) {
        invisible(dev.off())
    }
}

#' Plot mutational spectra
#'
#' \code{plot_spectrum} generates plots of one or more mutational spectra, which can be either
#' mutational catalogues or mutational signatures. If provided with multiple spectra, it produces
#' one plot per spectrum. If the spectra contain values greater than 1, the values will be
#' interpreted as mutation counts (as in a catalogue); otherwise, they will be interpreted as
#' mutation probabilities (as in a signature).
#' @param spectra This can be a numeric vector with one element per mutation type, a numeric matrix
#' with one row per signature/catalogue and one column per mutation type, or a list of signature
#' matrices as produced by \code{\link{retrieve_pars}}. In the latter case, HPD intervals will be 
#' included in the plots. Row names will be used as the catalogue/signature names.
#' @param pdf_path Character indicating the path to an optional output PDF file for the plots. The
#' PDF dimensions and graphical parameters are automatically set to appropriate values, but custom
#' dimensions can be specified via the arguments \code{pdf_width} and \code{pdf_height}.
#' @param pdf_width Numeric indicating the width of the output PDF, in inches (default is 24).
#' Only used if \code{pdf_path} is provided.
#' @param pdf_height Numeric indicating the height of the output PDF, in inches (default is 8).
#' Only used if \code{pdf_path} is provided.
#' @param name Character indicating a name to include in the plot title; useful when plotting a
#' single spectrum.
#' @param max_y Numeric indicating an optional upper limit for the vertical axis.
#' @param colors Character vector of custom color names or hexadecimal codes to use for the spectrum
#' bars. Only used if the number of mutation types in the spectrum is not 96 or 192. Must contain
#' either a single value, or as many values as the number of mutation types in the spectrum.
#' @param boxes Logical indicating whether boxes should be drawn around the plots (default is
#' \code{TRUE}).
#' @examples
#' \dontrun{
#' # Load example mutational catalogues
#' data("counts_21breast")
#'
#' # Plot catalogues
#' plot_spectrum(counts_21breast, pdf_path = "Catalogues.pdf")
#'
#' # Extract signatures using the Poisson model
#' samples <- extract_signatures(counts_21breast, nsignatures = 2, model = "poisson",
#'                               opportunities = "human-genome", iter = 800)
#'
#' # Retrieve extracted signatures
#' sigs <- retrieve_pars(samples, "signatures")
#'
#' # Plot signatures
#' plot_spectrum(sigs, pdf_path = "Signatures.pdf")
#' }
#' @importFrom "grDevices" cairo_pdf
#' @export
plot_spectrum <- function(spectra, pdf_path = NULL, pdf_width = 24, pdf_height = 8,
                          name = NULL, max_y = NULL, colors = NULL, boxes = TRUE) {
    # Fetch HPD interval values, if present
    if (is.list(spectra) & "mean" %in% names(spectra)) {
        spec <- to_matrix(spectra$mean)
        lwr <- to_matrix(spectra$lower)
        upr <- to_matrix(spectra$upper)
    }
    else {
        spec <- to_matrix(spectra)
        lwr <- NULL
        upr <- NULL
    }

    NCAT <- ncol(spec)      # number of categories
    NSAMP <- nrow(spec)     # number of samples
    strand <- NCAT == 192   # strand bias indicator (logical)
    counts <- any(spec > 1) # count data indicator

    # Plot each spectrum
    TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    STRANDCOL <- c("deepskyblue3", "red3")
    BACKCOL <- c("#00BFFF33", "#00000033", "#EE2C2C33", "#C2C2C24D", "#A2CD5A4D", "#EEB4B44D")
    LINECOL <- "gray60"
    XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    BACKLIM <- c(0, 46.8, 93.2, 139.55, 186, 232.35, 279.2)

    if (!is.null(pdf_path)) {
        cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, onefile = TRUE)
        if (ncol(spec) %in% c(96, 192)) {
            par(mar = c(5.5, 7, 7.5, 2))
        }
        else {
            par(mar = c(10, 7, 7.5, 2))
        }
    }

    # Generic spectrum (NCAT!={96,192})
    if (!(ncol(spec) %in% c(96, 192))) {
        if (is.null(colnames(spec))) {
            types <- paste("Mut. type", 1:ncol(spec))
        }
        else {
            types <- colnames(spec)
        }

        for (i in 1:NSAMP) {
            if (is.null(max_y)) {
                FACTOR <- 1.05
                samp_max_y <- ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR)
            }
            else {
                samp_max_y <- max_y
            }
            if (is.null(colors)) {
                colors = "orangered3"
            }
            else {
                if ((length(colors) > 1) & (length(colors) != ncol(spec))) {
                    stop("'colors' must contain either a single value, or one value per mutation type.")
                }
            }
            # Plot spectrum bars
            bars <- barplot(spec[i,], names.arg = types, col = colors, mgp = c(3, 0.8, 0),
                            border = "white", las = 2, cex.names = 1,
                            ylim = c(0, samp_max_y), yaxt = "n")
            # Plot axis
            if (counts) {
                axis(side = 2, cex.axis = 1.9, lwd = 2)
                label <- "Mutations"
                n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
            }
            else {
                axis(side = 2, cex.axis = 1.9, lwd = 2)
                label <- "Mutation probability"
                n_text <- ""
            }
            if (is.null(name)) {
                nme <- rownames(spec)[i]
            }
            else {
                nme <- name
            }
            mtext(label, side = 2, cex = 2.4, line = 3.5)
            title(paste0(nme, n_text), line = 4, cex.main = 2.5)
            # Plot HPD intervals
            if (!is.null(lwr)) {
                arrows(bars, upr[i,],
                       bars, lwr[i,],
                       length = 0, lwd = 3, col = LINECOL)
            }
            # Plot box
            if (boxes) {
                box(lwd = 2)
            }
        }
    }
    else {

        # Standard spectrum (NCAT=96)
        if (!strand) {
            for (i in 1:NSAMP) {
                if (is.null(max_y)) {
                    FACTOR <- 1.095
                    samp_max_y <- max(0.05,
                                      ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
                }
                else {
                    samp_max_y <- max_y
                }
                if (boxes) {
                    xlim <- c(-0.105, 115.5)
                }
                else {
                    xlim <- c(-1, 116)
                }
                # Plot spectrum bars
                bars <- barplot(spec[i, ],
                                names.arg = substr(mut_types(), 1, 3), mgp = c(3, 0.8, 0),
                                col = rep(COLORS, each = 16), border = "white",
                                las = 2, ylim = c(0, samp_max_y), xlim = xlim,
                                yaxt = "n", cex.names = 1.6, xaxs = "i", family = "mono")
                # Highlight trinucleotide middle bases
                for (j in 1:length(COLORS)) {
                    idx <- ((j-1) * 16 + 1):(j * 16)
                    axis(side = 1, at = bars[idx], tick = FALSE, cex.axis = 1.6,
                         mgp = c(3, 0.8, 0), las = 2, family = "mono", font = 2,
                         col.axis = COLORS[j], labels = paste0(" ", substr(mut_types()[idx], 2, 2), " "))
                }
                # Plot axis
                if (counts) {
                    axis(side = 2, cex.axis = 1.9, lwd = 2)
                    label <- "Mutations"
                    n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
                }
                else {
                    axis(side = 2, at = seq(0, samp_max_y, 0.05), cex.axis = 1.9, lwd = 2)
                    label <- "Mutation probability"
                    n_text <- ""
                }
                if (is.null(name)) {
                    nme <- rownames(spec)[i]
                }
                else {
                    nme <- name
                }
                mtext(label, side = 2, cex = 2.4, line = 3.5)
                title(paste0(nme, n_text), line = 4, cex.main = 2.5)
                # Plot HPD intervals
                if (!is.null(lwr)) {
                    arrows(bars, upr[i,],
                           bars, lwr[i,],
                           length = 0, lwd = 3, col = LINECOL)
                }
                # Plot mutation type labels
                text(x = (XL + XR) / 2, y = 1.055 * samp_max_y,
                     labels = TYPES, cex = 2.4, xpd = TRUE)
                rect(xleft = XL, xright = XR, ybottom = 0.945 * samp_max_y,
                     ytop = samp_max_y, col = COLORS, border = "white")
                # Plot box
                if (boxes) {
                    box(lwd = 2)
                }
            }
        }

        # Strand-wise spectrum (NCAT=192)
        else {
            for (i in 1:NSAMP) {
                if (is.null(max_y)) {
                    FACTOR <- 1.095
                    samp_max_y <- max(0.05,
                                      ifelse(is.null(upr), max(spec[i,]) * FACTOR, max(upr[i,]) * FACTOR))
                }
                else {
                    samp_max_y <- max_y
                }
                if (boxes) {
                    xlim <- c(0, 279.2)
                }
                else {
                    xlim <- c(-3, 280)
                }
                # Plot background panes and mutation type labels
                barplot(rbind(spec[i, 1:(NCAT/2)], spec[i, (NCAT/2+1):NCAT]), beside = TRUE,
                        col = NA, border = NA, space = c(0.1, 0.8), xaxs = "i", yaxt = "n",
                        xaxt = "n", ylim = c(0, samp_max_y), xlim = xlim)
                for (j in 1:length(COLORS)) {
                    rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0,
                         ytop = samp_max_y, col = BACKCOL[j], border = "white")
                    text(x = (BACKLIM[j] + BACKLIM[j+1]) / 2, y = 1.055 * samp_max_y,
                         labels = TYPES[j], cex = 2.4, xpd = TRUE)
                    rect(xleft = BACKLIM[j], xright = BACKLIM[j+1], ybottom = 0.945 * samp_max_y,
                         ytop = samp_max_y, col = COLORS[j], border = "white")
                }
                # Plot legend
                legend("topright", bty = "n", inset = c(0.016, 0.03), 
                       legend = c("Transcribed", "Untranscribed"), 
                       cex = 2.1, fill = NA, border = NA)
                legend("topright", bty = "n", inset = c(0.115, 0.03), pch=15, pt.cex=3.75,
                       col=STRANDCOL, legend = c("", ""), cex = 2.1)
                # Plot spectrum bars
                bars <- barplot(rbind(spec[i, 1:(NCAT/2)],
                                      spec[i, (NCAT/2+1):NCAT]),
                                names.arg = substr(mut_types(), 1, 3), beside = TRUE,
                                space = c(0.1, 0.8), mgp = c(3, 0.8, 0), las = 2,
                                col = STRANDCOL, border = "white", yaxt = "n",
                                cex.names = 1.6, xaxs = "i", family = "mono", add = TRUE)
                # Highlight trinucleotide middle bases
                for (j in 1:length(COLORS)) {
                    idx <- ((j-1) * 16 + 1):(j * 16)
                    axis(side = 1, tick = FALSE, at = colMeans(bars[, idx]),
                         cex.axis = 1.6, mgp = c(3, 0.8, 0), las = 2,
                         family = "mono", font = 2, col.axis = COLORS[j],
                         labels = paste0(" ", substr(mut_types()[idx], 2, 2), " "))
                }
                # Plot axis
                if (counts) {
                    axis(side = 2, cex.axis = 1.9, lwd = 2)
                    label <- "Mutations"
                    n_text <- paste0(" (", prettyNum(sum(spec[i,]), big.mark = ","), " mutations)")
                }
                else {
                    axis(side = 2, at = seq(0, samp_max_y, 0.05), cex.axis = 1.9, lwd = 2)
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
                mtext(label, side = 2, cex = 2.4, line = 3.5)
                title(paste0(nme, n_text), line = 4, cex.main = 2.5)
                # Plot HPD intervals
                if (!is.null(lwr)) {
                    bars <- as.numeric(t(bars))
                    arrows(bars, upr[i,],
                           bars, lwr[i,],
                           length = 0, lwd = 2.5, col = LINECOL)
                }
                # Plot box
                if (boxes) {
                    box(lwd = 2)
                }
            }
        }
    }

    if (!is.null(pdf_path)) {
        invisible(dev.off())
    }
}

#' Plot signature exposures
#'
#' \code{plot_exposures} plots the distribution of signature exposures across samples.
#' @param mcmc_samples List with two elements named \code{`data`} and \code{`results`}, produced via
#' \code{\link{fit_signatures}}, \code{\link{extract_signatures}}, or
#' \code{\link{fit_extract_signatures}}. This is the preferred option for supplying data and
#' results, but can be replaced by the combination of arguments \code{counts}, \code{exposures} and
#' \code{signature_names}.
#' @param pdf_path Character indicating the path to an optional output PDF file for the plots. The
#' PDF dimensions and graphical parameters are automatically set to appropriate values, but custom
#' dimensions can be specified via the arguments \code{pdf_width} and \code{pdf_height}.
#' @param counts Numeric matrix of observed mutation counts, with one row per sample and
#' one column per mutation type. Only needed if \code{mcmc_samples} is not provided.
#' @param exposures Either a numeric matrix of signature exposures, with one row per sample and one
#' column per signature, or a list of matrices generated via \code{\link{retrieve_pars}}. Only
#' needed if \code{mcmc_samples} is not provided.
#' @param signature_names Character vector containing the name of each signature. Only needed if
#' \code{mcmc_samples} is not provided and the exposures were obtained via signature fitting
#' (rather than extraction).
#' @param thresh Numeric indicating the minimum threshold for the lower HPD limits of signature
#' exposures (default is 0.01). Exposures with a lower HPD limit below this value will be shown in
#' grey.
#' @param hpd_prob Numeric value in the interval (0, 1), indicating the desired probability content
#' of HPD intervals (default is 0.95).
#' @param pdf_width Numeric indicating the width of the output PDF, in inches (default is 24).
#' Only used if \code{pdf_path} is provided.
#' @param pdf_height Numeric indicating the height of the output PDF, in inches (default is 10).
#' Only used if \code{pdf_path} is provided.
#' @param margin_bottom Numeric indicating the width of the bottom margin, in inches (default
#' is 10.5).
#' @param legend_pos Character indicating the position of the legend. Admits values \code{"top"},
#' \code{"bottom"}, \code{"center"}, \code{"left"}, \code{"right"}, \code{"topleft"},
#' \code{"topright"}, \code{"bottomleft"} and \code{"bottomright"} (default is \code{"topleft"}).
#' @param legend_cex Numeric indicating the relative size of the legend (default is 2).
#' @param cex_names Numeric indicating the relative size of sample labels (default is 1.9).
#' @param sig_color_palette Character vector of custom color names or hexadecimal codes to use for
#' each signature in exposure and reconstruction plots. Must have at least as many elements as the
#' number of signatures.
#' @importFrom "graphics" arrows axis barplot legend lines mtext par plot points rect text title
#' @importFrom "grDevices" pdf dev.off rgb
#' @examples
#' \dontrun{
#' # Load example mutational catalogues and COSMIC signatures
#' data("counts_21breast")
#' data("cosmic_signatures_v2")
#'
#' # Fit signatures and retrieve exposures
#' samples <- fit_signatures(counts_21breast, cosmic_signatures_v2)
#' exposures <- retrieve_pars(samples, "exposures")
#'
#' # Plot exposures using MCMC samples
#' plot_exposures(mcmc_samples = samples, pdf_path = "Exposures.pdf")
#'
#' # Plot exposures using retrieved exposures matrix
#' plot_exposures(counts = counts_21breast, exposures = exposures,
#'                signature_names = rownames(cosmic_signatures_v2),
#'                pdf_path = "Exposures.pdf")
#' }
#' @importFrom "grDevices" cairo_pdf
#' @importFrom "graphics" segments
#' @export
plot_exposures <- function(mcmc_samples = NULL, pdf_path = NULL, counts = NULL, exposures = NULL,
                           signature_names = NULL, thresh = 0.01, hpd_prob = 0.95, pdf_width = 24,
                           pdf_height = 10, margin_bottom = 10.5, legend_pos = "topleft",
                           legend_cex = 2, cex_names = 1.9, sig_color_palette = NULL) {
    if (is.null(mcmc_samples) & (is.null(counts) | is.null(exposures))) {
        stop("Either 'mcmc_samples', or both 'counts' and 'exposures', must be provided.")
    }
    if (!is.null(mcmc_samples)) {
        counts <- mcmc_samples$data$counts_real
        exposures <- retrieve_pars(mcmc_samples, "exposures", hpd_prob = hpd_prob)
        lwr <- to_matrix(exposures$lower)
        upr <- to_matrix(exposures$upper)
    }
    else if (is.list(exposures) & "mean" %in% names(exposures)) {
        lwr <- to_matrix(exposures$lower)
        upr <- to_matrix(exposures$upper)
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
        if (pdf_width == 24) {
            pdf_width <- max(pdf_width, NSAMP * 0.13)
        }
        cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, onefile = TRUE)
        par(mar = c(margin_bottom, 7.5, 7.5, 0))
    }

    # If >1 sample: plot average exposures
    if (NSAMP > 1) {
        exposures_global <- colMeans(exposures)
        if (!is.null(lwr)) {
            lwr_global <- colMeans(lwr)
            upr_global <- colMeans(upr)
            max_y <- max(upr_global)
        }
        else {
            max_y <- max(exposures_global)
        }
        colours <- rep("skyblue3", NSIG)
        if (!is.null(lwr)) {
            colours[lwr_global < thresh] <- "grey"
        }
        bars <- barplot(exposures_global, col = colours, border = NA,
                        cex.names = 1e-20, cex.main = 2.3, ylim = c(0, max_y), axes = F,
                        main = "Mean signature exposures across sample set")
        text(x = bars, y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
             labels = names(exposures_global), cex = cex_names, srt = 45, adj = 1, xpd = TRUE)
        axis(side = 2, cex.axis = 1.9, lwd = 2, line = -2.5, las = 2)
        mtext("Mutation fraction", side = 2, cex = 2.4, line = 3)
        if (!is.null(lwr)) {
            arrows(bars, upr_global,
                   bars, lwr_global,
                   length = 0, lwd = 4, col = "gray50")
        }
    }

    # Plot exposures for each sample
    if (!is.null(upr)) {
        max_y <- max(upr)
    }
    else {
        max_y <- max(exposures)
    }
    for (i in 1:NSAMP) {
        colours <- rep("skyblue3", NSIG)
        if (!is.null(lwr)) {
            colours[lwr[i, ] < thresh] <- "grey"
        }
        bars <- barplot(exposures[i, ], col = colours, border = NA,
                        cex.names = 1e-20, cex.main = 2.3, ylim = c(0, max_y), axes = F,
                        main = paste("Signature exposures in", rownames(exposures)[i]))
        text(x = bars, y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
             labels = colnames(exposures), cex = cex_names, srt = 45, adj = 1, xpd = TRUE)
        axis(side = 2, cex.axis = 1.9, lwd = 2, line = -2.5, las = 2)
        mtext("Mutation fraction", side = 2, cex = 2.4, line = 3)
        if (!is.null(lwr)) {
            segments(x0 = bars, y0 = upr[i, ], y1 = lwr[i, ], lwd = 4, col = "gray50")
        }
    }

    # If >1 sample: plot exposures across samples
    if (NSAMP > 1) {
        # Plot absolute exposures
        muts <- rowSums(counts)
        exposures_abs <- exposures * muts
        bars <- barplot(t(exposures_abs), col = sigcols,
                        cex.names = 1e-20, cex.main = 2.3, axes = FALSE,
                        main = "Signature exposures per sample (absolute)")
        text(x = bars, y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
             labels = rownames(exposures), cex = cex_names, srt = 45, adj = 1, xpd = TRUE)
        axis(side = 2, cex.axis = 1.9, lwd = 2, line = -2.5)
        mtext("Mutations", side = 2, cex = 2.4, line = 3)

        # Legend
        # Expand legend box horizontally if there are many signatures
        LEGENDCOLS <- max(2, ceiling(NSIG / 10))
        legend(legend_pos, bty = "n", ncol = LEGENDCOLS, xpd = TRUE, inset = c(0.03, -0.04),
               fill = sigcols, border = sigcols, legend = colnames(exposures), cex = legend_cex)

        # Plot relative exposures
        bars <- barplot(t(exposures), col = sigcols,
                        cex.names = 1e-20, cex.main = 2.3, axes = FALSE,
                        main = "Signature exposures per sample (relative)")
        text(x = bars, y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
             labels = rownames(exposures), cex = cex_names, srt = 45, adj = 1, xpd = TRUE)
        axis(side = 2, cex.axis = 1.9, lwd = 2, line = -2.5, las = 2)
        mtext("Mutation fraction", side = 2, cex = 2.4, line = 3)
    }

    if (!is.null(pdf_path)) {
        invisible(dev.off())
    }
}
