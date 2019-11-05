## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library(sigfit))
par(mar = c(6, 4, 6, 4))

## ----devtools_instructions, eval=FALSE-----------------------------------
#  devtools::install_github("kgori/sigfit", build_opts = c("--no-resave-data", "--no-manual"))

## ----fetch---------------------------------------------------------------
library(sigfit)
data("cosmic_signatures_v2")

## ----sim-----------------------------------------------------------------
set.seed(1)
probs <- c(0.4, 0.3, 0.2, 0.1) %*% as.matrix(cosmic_signatures_v2[c(1, 3, 7, 11), ])
mutations <- matrix(rmultinom(1, 20000, probs), nrow = 1)
colnames(mutations) <- colnames(cosmic_signatures_v2)

## ----plotsim, fig.width=20, fig.height=7, out.width="100%", echo=-1------
par(mar = c(4.5,5.5,6.5,1))
plot_spectrum(mutations, name = "Simulated counts")

## ----fitting, warning=FALSE----------------------------------------------
mcmc_samples_fit <- fit_signatures(counts = mutations, 
                                   signatures = cosmic_signatures_v2,
                                   iter = 2000, 
                                   warmup = 1000, 
                                   chains = 1, 
                                   seed = 1756)

## ----retrieve_exp--------------------------------------------------------
exposures <- retrieve_pars(mcmc_samples_fit, 
                           par = "exposures", 
                           hpd_prob = 0.90)
names(exposures)
exposures$mean

## ----plot_exp, fig.width=17, fig.height=7, out.width='100%', fig.align="center", echo=-1----
par(mar=c(8,5,3.5,0))
plot_exposures(mcmc_samples = mcmc_samples_fit)

## ----refitting, eval=FALSE-----------------------------------------------
#  selected_signatures <- c(1, 3, 7, 11)
#  mcmc_samples_fit_2 <- fit_signatures(mutations,
#                                       cosmic_signatures_v2[selected_signatures, ],
#                                       iter = 2000,
#                                       warmup = 1000,
#                                       chains = 1,
#                                       seed = 1756)

## ----reconstruct, fig.width=25, fig.height=18.5, out.width='100%', warning=FALSE, results="hide", echo=-1----
par(mar=c(5,6,6.5,1))
plot_reconstruction(mcmc_samples = mcmc_samples_fit)

## ----plot_all, eval=FALSE------------------------------------------------
#  plot_all(mcmc_samples = mcmc_samples_fit,
#           out_path = "your/output/dir/here",
#           prefix = "Fitting")

## ----load_mutations------------------------------------------------------
library(sigfit)
data("variants_21breast")
head(variants_21breast)

## ----show_samples--------------------------------------------------------
unique(variants_21breast[, 1])

## ----build_catalogues----------------------------------------------------
counts_21breast <- build_catalogues(variants_21breast)
dim(counts_21breast)
counts_21breast[1:5, 1:8]

## ----plot_spectra, fig.width=22, fig.height=23, out.width='100%', fig.align="center", echo=-1----
par(mar = c(5,6,7,2))
par(mfrow = c(7, 3))
plot_spectrum(counts_21breast)

## ----extraction, eval=FALSE----------------------------------------------
#  mcmc_samples_extr <- extract_signatures(counts = counts_21breast,
#                                          nsignatures = 2:7,
#                                          iter = 1000,
#                                          seed = 1756)

## ----plot_gof_silent, echo=FALSE, fig.width=9, fig.height=6, out.width="100%"----
## Plot precalculated GOF in order to avoid running the model
data("sigfit_vignette_data", package = "sigfit")
plot(nS, gof, type = "o", lty = 3, pch = 16, col = "dodgerblue4",
     main = paste0("Goodness-of-fit (", stat, ")\nmodel: multinomial"),
     xlab = "Number of signatures", 
     ylab = paste0("Goodness-of-fit (", stat, ")"))
points(nS[best], gof[best], pch = 16, col = "orangered", cex = 1.1)
cat("Estimated best number of signatures:", nS[best], "\n")

## ----extr_names, eval=FALSE----------------------------------------------
#  names(mcmc_samples_extr)

## ----extr_names_silent, echo=FALSE---------------------------------------
print(c("nsignatures=1", "nsignatures=2", "nsignatures=3", "nsignatures=4", "nsignatures=5", "nsignatures=6", "nsignatures=7", "best"))
signatures <- extr_signatures

## ----retrieve_sigs, eval=FALSE-------------------------------------------
#  ## Note: mcmc_samples_extr[[N]] contains the extraction results for N signatures
#  signatures <- retrieve_pars(mcmc_samples_extr[[4]],
#                              par = "signatures")

## ----show_signames-------------------------------------------------------
rownames(signatures$mean)

## ----plot_sigs, warning=FALSE, fig.width=25, fig.height=10, out.width='100%', fig.align="center", echo=-1----
par(mar = c(6,7,6,1))
par(mfrow = c(2, 2))
plot_spectrum(signatures)

## ----match_sigs----------------------------------------------------------
data("cosmic_signatures_v2")
match_signatures(signatures, cosmic_signatures_v2)

## ----refitting2, eval=FALSE----------------------------------------------
#  signatures <- retrieve_pars(mcmc_samples_extr_2, "signatures")
#  mcmc_samples_refit <- fit_signatures(counts = counts_21breast,
#                                       signatures = signatures,
#                                       iter = 2000,
#                                       warmup = 1000)
#  exposures <- retrieve_pars(mcmc_samples_refit, "exposures")

## ----strandwise, fig.width=22, fig.height=7, out.width="100%", echo=-1----
par(mar = c(4.5,5.5,6.5,1))
# Load and plot a strand-wise catalogue
data("counts_88liver_strand")
plot_spectrum(counts_88liver_strand[1, ], name = rownames(counts_88liver_strand)[1])

## ----generic, fig.width=22, fig.height=8, out.width="100%", echo=-1------
par(mar = c(6.5,5.5,6.5,1)); set.seed(0xC0FFEE)
# Create and plot an arbitrary catalogue
counts <- as.numeric(rmultinom(1, 5000, runif(100)))
plot_spectrum(counts, name = "Arbitrary catalogue")

## ----convert_sigs, eval=F------------------------------------------------
#  # (Following from the signature extraction example presented above)
#  # Retrieve genome-derived signatures
#  genome_signatures <- retrieve_pars(mcmc_samples_extr[[4]],
#                                     par = "signatures")
#  
#  # Apply exome mutational opportunities
#  exome_signatures <- convert_signatures(genome_signatures,
#                                         opportunities_from = "human-genome",
#                                         opportunities_to = "human-exome")
#  
#  par(mfrow = c(2, 1))
#  plot_spectrum(genome_signatures$mean[4,],
#                name = "Signature D, Genome-relative probabilities")
#  plot_spectrum(exome_signatures[4,],
#                name = "Signature D, Exome-relative probabilities")

## ----convert_sigs_silent, echo=F, fig.width=22, fig.height=14, out.width="100%"----
exome_signatures <- convert_signatures(signatures, opportunities_from = "human-genome", opportunities_to = "human-exome")
par(mfrow = c(2, 1), mar = c(5,5.5,6.5,1))
plot_spectrum(signatures$mean[4,], name="Signature D, Genome-relative probabilities")
plot_spectrum(exome_signatures[4,], name="Signature D, Exome-relative probabilities")

## ----multichain, eval=F--------------------------------------------------
#  # Load mutation counts
#  data("counts_21breast")
#  
#  # Set number of cores for multi-chain sampling
#  ncores <- parallel::detectCores()
#  options(mc.cores = ncores)
#  
#  # Run a short single chain to find initial parameter values
#  ext_init <- extract_signatures(counts_21breast,
#                                 nsignatures = 4,
#                                 iter = 1000,
#                                 chains = 1)
#  
#  # Obtain a list of initial parameter values
#  inits <- get_initializer_list(ext_init, chains = ncores)
#  
#  # Run multiple chains in parallel, initialised at the values found in the short run
#  # (make sure to use the same number of chains as in `get_initializer_list`)
#  ext_multichain <- extract_signatures(counts_21breast,
#                                       nsignatures = 4,
#                                       iter = 5000,
#                                       chains = ncores,
#                                       init = inits)

## ------------------------------------------------------------------------
    data("cosmic_signatures_v2")
    data("cosmic_signatures_v3")
    data("cosmic_signatures_v3_strand")

## ------------------------------------------------------------------------
    data("variants_21breast")
    data("counts_21breast")

## ------------------------------------------------------------------------
    data("variants_119breast")
    data("variants_119breast_strand")
    data("counts_119breast")
    data("counts_119breast_strand")

## ------------------------------------------------------------------------
    data("variants_88liver")
    data("variants_88liver_strand")
    data("counts_88liver")
    data("counts_88liver_strand")

## ------------------------------------------------------------------------
    data("methylation_50breast")
    data("methylation_27normal")

