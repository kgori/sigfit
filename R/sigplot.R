### Helper functions (not exported)

#' Extract signature data from extracted rstan samples array
#' @param sig item such as stanfit_extract$signatures[, 1, ]
#' @importFrom "coda" HPDinterval as.mcmc
signature_dataframe <- function(sig) {
    n <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A",
           "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C",
           "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G",
           "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
           "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A",
           "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C",
           "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G",
           "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
           "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A",
           "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C",
           "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
    hpd <- HPDinterval(as.mcmc(sig))
    data.frame(context = factor(n, levels = n),
               lower = hpd[, "lower"],
               rate = colMeans(sig),
               upper = hpd[, "upper"],
               mutation = factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), times = rep(16, 6))))
}

#' Extract exposures data from extracted rstan samples array
#' @param exp item such as stanfit_extract$signatures[, 1, ]
#' @importFrom "coda" HPDinterval as.mcmc
exposure_dataframe <- function(exp) {
    hpd <- HPDinterval(as.mcmc(exp))
    data.frame(signature = 1:ncol(exp),
               lower = hpd[, "lower"],
               exposure = colMeans(exp),
               upper = hpd[, "upper"])
}


### Exported functions

#' Returns a list of signature spectrum plots, with error bars.
#' @param samples Sampling object such as returned by extract_signatures.
#' @examples
#' Generate a list of 3 plots for 3 extracted signatures
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "nmf") 
#' plots <- plot_signatures(samples)
#' m <- marrangeGrob(plots, nrow = 1, ncol = 3) # use gridExtra to plot all on same page
#' @importFrom "rstan" extract
#' @importFrom "ggplot2" ggplot geom_col geom_errorbar scale_fill_manual theme theme_bw element_text ggtitle ylim
#' @export
ggplot_signatures <- function(samples) {
    signatures <- extract(samples)$signatures
    plots <- list()
    max.y <- 0
    for (i in 1:dim(signatures)[2]) {
        df <- signature_dataframe(signatures[, i,])
        max.y <- max(max.y, max(df$upper))
        colours <- c("deepskyblue", "black", "firebrick2",
                     "gray76", "darkolivegreen3", "rosybrown2")
        plots[[i]] <- ggplot(df, aes(context, rate)) +
            geom_col(aes(fill = mutation)) +
            geom_errorbar(
                aes(ymin = lower, ymax = upper),
                colour = "grey40",
                alpha = 0.7,
                width = 0.5
            ) +
            scale_fill_manual(values = colours) +
            theme_bw() +
            theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0.5,
                size = 8
            )) +
            ggtitle(paste("Extracted signature", i))
    }
    for (i in 1:length(plots)) {
        plots[[i]] <- plots[[i]] + ylim(0, max.y)
    }
    plots
}


#' Returns a list of normalised exposure plots, with error bars.
#' @param samples Sampling object such as returned by extract_signatures.
#' @examples
#' Generate a list of exposure plots for all samples, showing 3 extracted signatures
#' samples <- extract_signatures(mycounts, nsignatures = 3, method = "nmf") 
#' plots <- plot_exposures(samples)
#' m <- marrangeGrob(plots, nrow = 1, ncol = 3) # use gridExtra to plot all on same page
#' @importFrom "rstan" extract
#' @importFrom "ggplot2" ggplot geom_col geom_errorbar scale_fill_manual guides theme_bw ggtitle ylim
#' @export
ggplot_exposures <- function(samples) {
    exposures <- extract(samples)$exposures
    colours <- c("dodgerblue", "grey80")
    max.y <- 0
    for (i in 1:dim(exposures)[2]) {
        df <- exposure_dataframe(exposures[, i,])
        max.y <- max(max.y, df$upper)
        plots[[i]] <- ggplot(df, aes(signature, exposure)) +
            geom_col(aes(fill = lower < 0.001), colour = "black") +
            geom_errorbar(aes(ymin = lower, ymax = upper),
                          colour = "black",
                          width = 0.1) +
            scale_fill_manual(values = colours) +
            theme_bw() + 
            guides(fill = FALSE) +
            ggtitle(paste0("Sample ", i, ": Relative signature exposures"))
    }
    for (i in 1:dim(exposures)[2]) {
        plots[[i]] <- plots[[i]] + ylim(0, max.y)
    }
    plots
}
