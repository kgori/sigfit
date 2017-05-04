## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library(sigfit))

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("kgori/sigfit")

## ----fetch---------------------------------------------------------------
sigs <- sigfit::fetch_cosmic_data()

## ----sim-----------------------------------------------------------------
set.seed(1)
probs <- as.matrix(sigs[, c(1, 5, 7, 11)]) %*% c(0.4, 0.3, 0.2, 0.1)
mutations <- as.vector(rmultinom(1, 20000, probs))
names(mutations) <- rownames(sigs)

## ----fig.width=7.5, fig.height=5, plotsim--------------------------------
barplot(mutations, cex.names = 0.4, las = 2)

## ----include=FALSE-------------------------------------------------------
set.seed(1)
mcmc_samples <- sigfit::run_sampling(counts = mutations, signatures = sigs,
                                       iter = 1500, warmup = 500, seed = 1)

## ----eval=FALSE----------------------------------------------------------
#  set.seed(1)
#  mcmc_samples <- sigfit::run_sampling(counts = mutations, signatures = sigs,
#                                         iter = 1500, warmup = 500, seed = 1)

## ----fig.width=7.5, fig.height=5, echo=FALSE-----------------------------
anon <- suppressWarnings(sigfit::plot_spectrum(mcmc_samples, names.arg = names(mutations),
                                                 cex.names = 0.4, las = 2))

## ----eval=FALSE----------------------------------------------------------
#  sigfit::plot_spectrum(mcmc_samples, names.arg = names(mutations),
#                          cex.names = 0.4, las = 2)

## ----fig.width=7.5, fig.height=5, exposures------------------------------
sigfit::plot_exposures(mcmc_samples, prob = 0.9, thresh = 0.01,
                         cex.names = 0.7, las = 2)

