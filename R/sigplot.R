#' Plots the mutational signatures in a set of samples
#' @importFrom "rstan" extract
#' @export
plot_signatures <- function(samples) {
    sigs <- extract(samples, pars = "signatures")
    meansigs <- apply(sigs[[1]], c(2, 3), mean)
    
    for (i in 1:nrow(meansigs)) {
        sigplot(meansigs[i,], titletext = paste("Signature in sampled position", i))
    }
}
    
#' Plots a mutational spectrum
#' @param spectrum Vector of mutation rates
#' @param titletext A string to use as the title of the plot
sigplot <- function(spectrum, titletext = "") {
    colors <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    mut <- c('C>A','C>G','C>T','T>A','T>C','T>G')
    tickfreq <- 0.05 # Plot a tick on the y-axis this frequently
    # Observed spectrum
    names(spectrum) <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A",
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
    height <- max(max(spectrum), 0.25)
    # cat(height)
    barplot(spectrum, 
            col = rep(colors, each = 16), 
            border = "white", 
            family = "mono",
            yaxt = "n", 
            ylim = c(0, height), 
            xlim = c(-1, 116), 
            xaxs = "i",
            cex = .6,  # <- x axis label size
            cex.axis = 1.5, 
            cex.lab = 1.7, 
            las = 2)
    axis_max <- ceiling(height/tickfreq)*tickfreq
    axis(side = 2, at = seq(0,axis_max,tickfreq), las = 2, cex.axis = 1.25)
    mtext("Mutation probability", side = 2, cex = 1.7, line = 4.5)
    title(titletext, line = 2, cex.main = 1)
    xleft = c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
    xright = c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
    rect(xleft = xleft, xright = xright, ybottom = axis_max, ytop = axis_max+0.01, 
         col = colors, border = "white", lty = par("lty"), lwd = par("lwd"))
    text(x = (xleft+xright)/2, y = axis_max-0.01, labels = mut, cex = 1.25)
}
