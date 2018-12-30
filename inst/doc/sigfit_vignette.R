## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library(sigfit))
par(mar = c(6, 4, 6, 4))

## ----devtools_instructions, eval=FALSE-----------------------------------
#  devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)

## ----fetch---------------------------------------------------------------
data("cosmic_signatures", package = "sigfit")

## ----sim-----------------------------------------------------------------
set.seed(1)
probs <- c(0.4, 0.3, 0.2, 0.1) %*% as.matrix(cosmic_signatures[c(1, 3, 7, 11), ])
mutations <- matrix(rmultinom(1, 20000, probs), nrow = 1)
colnames(mutations) <- colnames(cosmic_signatures)

## ----plotsim, fig.width=17, fig.height=7, out.width="100%", echo=-1------
par(mar = c(4.5,5.5,6.5,1))
sigfit::plot_spectrum(mutations, name = "Simulated counts")

## ----fitting, warning=FALSE----------------------------------------------
mcmc_samples_fit <- sigfit::fit_signatures(counts = mutations, 
                                           signatures = cosmic_signatures,
                                           iter = 2000, 
                                           warmup = 1000, 
                                           chains = 1, 
                                           seed = 1)

## ----retrieve_exp--------------------------------------------------------
exposures <- sigfit::retrieve_pars(mcmc_samples_fit, 
                                   par = "exposures", 
                                   hpd_prob = 0.90)
names(exposures)
exposures$mean

## ----plot_exp, fig.width=17, fig.height=7, out.width='100%', fig.align="center", echo=-1----
par(mar=c(8,5,3.5,0))
sigfit::plot_exposures(mcmc_samples = mcmc_samples_fit)

## ----reconstruct, fig.width=25, fig.height=18.5, out.width='100%', warning=FALSE, results="hide", echo=-1----
par(mar=c(5,6,6.5,1))
sigfit::plot_reconstruction(mcmc_samples = mcmc_samples_fit,
                            pdf_path = NULL)

## ----plot_all, eval=FALSE------------------------------------------------------------------
#  sigfit::plot_all(mcmc_samples = mcmc_samples_fit,
#                   out_path = "your/output/dir/here",
#                   prefix = "Fitting")

## ----load_mutations------------------------------------------------------------------------
data("variants_21breast", package = "sigfit")
head(variants_21breast)

## ----show_samples--------------------------------------------------------------------------
unique(variants_21breast[, 1])

## ----build_catalogues--------------------------------------------------------------
counts_21breast <- sigfit::build_catalogues(variants_21breast)
dim(counts_21breast)
counts_21breast[1:5, 1:8]

## ----plot_spectra, fig.width=22, fig.height=23, out.width='100%', fig.align="center", echo=-1----
par(mar = c(5,6,7,2))
par(mfrow = c(7, 3))
sigfit::plot_spectrum(counts_21breast)

## ----extraction, eval=FALSE--------------------------------------------------------
#  mcmc_samples_extr <- sigfit::extract_signatures(counts = counts_21breast,
#                                                  nsignatures = 2:7,
#                                                  iter = 1000,
#                                                  seed = 1)

## ----plot_gof_silent, echo=FALSE, fig.width=9, fig.height=6, out.width="100%"------
## Plot precalculated GOF in order to avoid running the model
data("sigfit_vignette_data", package = "sigfit")
plot(nS, gof, type = "o", lty = 3, pch = 16, col = "dodgerblue4",
     main = paste0("Goodness of fit (", stat, ")\nmodel: NMF"),
     xlab = "Number of signatures", 
     ylab = paste0("Goodness of fit (", stat, ")"))
points(nS[best], gof[best], pch = 16, col = "orangered", cex = 1.1)
cat("Estimated best number of signatures:", nS[best], "\n")

## ----extr_names, eval=FALSE--------------------------------------------------------
#  names(mcmc_samples_extr)

## ----extr_names_silent, echo=FALSE-------------------------------------------------
print(c("nsignatures=1", "nsignatures=2", "nsignatures=3", "nsignatures=4", "nsignatures=5", "nsignatures=6", "nsignatures=7", "best"))

## ----retrieve_sigs, eval=FALSE-----------------------------------------------------
#  ## Note: mcmc_samples_extr[[N]] contains the extraction results for N signatures
#  extr_signatures <- sigfit::retrieve_pars(mcmc_samples_extr[[4]],
#                                           par = "signatures")

## ----show_signames-----------------------------------------------------------------
rownames(extr_signatures$mean)

## ----plot_sigs, warning=FALSE, fig.width=25, fig.height=10, out.width='100%', fig.align="center", echo=-1----
par(mar = c(6,7,6,1))
par(mfrow = c(2, 2))
sigfit::plot_spectrum(extr_signatures)

## ----convert_sigs, eval=F----------------------------------------------------------
#  # Following from the signature extraction example above
#  genome_signatures <- sigfit::retrieve_pars(mcmc_samples_extr[[4]],
#                                             par = "signatures")
#  normalised_signatures <- sigfit::convert_signatures(genome_signatures,
#                                                      ref_opportunities = "human-genome",
#                                                      model_to = "emu")
#  exome_signatures <- sigfit::convert_signatures(normalised_signatures,
#                                                 ref_opportunities = "human-exome",
#                                                 model_to = "nmf")
#  
#  par(mfrow = c(2, 1))
#  sigfit::plot_spectrum(genome_signatures$mean[4,], name = "Signature D, Genome-relative probabilities")
#  sigfit::plot_spectrum(exome_signatures[4,], name = "Signature D, Exome-relative probabilities")

## ----convert_sigs_silent, echo=F, fig.width=20, fig.height=14, out.width="100%"----
genome_signatures <- extr_signatures
normalised_signatures <- sigfit::convert_signatures(genome_signatures, ref_opportunities="human-genome", model_to="emu")
exome_signatures <- sigfit::convert_signatures(normalised_signatures, ref_opportunities="human-exome", model_to="nmf")
par(mfrow = c(2, 1), mar = c(5,5.5,6.5,1))
sigfit::plot_spectrum(genome_signatures$mean[4,], name="Signature D, Genome-relative probabilities")
sigfit::plot_spectrum(exome_signatures[4,], name="Signature D, Exome-relative probabilities")

