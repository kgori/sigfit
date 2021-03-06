% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigfit_estimation.R
\name{fit_signatures}
\alias{fit_signatures}
\title{Fit mutational signatures}
\usage{
fit_signatures(
  counts,
  signatures,
  exp_prior = NULL,
  model = "multinomial",
  opportunities = NULL,
  ...
)
}
\arguments{
\item{counts}{Numeric matrix of observed mutation counts, with one row per sample and
one column per mutation type.}

\item{signatures}{Mutational signatures to be fitted; either a numeric matrix with one row per
signature and one column per mutation type, or a list of matrices generated via
\code{\link{retrieve_pars}}.}

\item{exp_prior}{Numeric matrix with one row per sample and one column per signature, to be used
as the Dirichlet priors for the signature exposures. Default priors are uniform.}

\item{model}{Name of the model to sample from. Admits character values \code{"multinomial"}
(default), \code{"poisson"}, \code{"negbin"}, \code{"normal"}, \code{"nmf"} (an alias for
\code{"multinomial"}), and \code{"emu"} (an alias for \code{"poisson"}).}

\item{opportunities}{Numeric matrix of optional mutational opportunities, with one row per sample
and one column per mutation type. It also admits character values \code{"human-genome"} or
\code{"human-exome"}, in which case the mutational opportunities of the reference human
genome/exome will be used for every sample.}

\item{...}{Additional arguments to be passed to \code{\link{rstan::sampling}}.}
}
\value{
A list with two elements:
\itemize{
 \item{\code{`data`}: list containing the input data supplied to the model.}
 \item{\code{`result`}: object of class stanfit, containing the output MCMC samples,
 as well as information about the model and the sampling process.}}
The model parameters (such as signature exposures) can be extracted from this
object using \code{\link{retrieve_pars}}.
}
\description{
\code{fit_signatures} performs MCMC sampling to fit a set of mutational signatures to a
collection of mutational catalogues and estimate the exposure of each catalogue to each signature.
Four models of signatures are available: multinomial, Poisson, normal and negative binomial. The
normal model can be used when \code{counts} contains continuous (non-integer) values, while the
negative binomial model is a more noise-robust version of the Poisson model.
}
\examples{
\dontrun{
# Load example mutational catalogues and COSMIC signatures
data("counts_21breast")
data("cosmic_signatures_v2")

# Fit signatures 1 to 4, using a custom prior that favors signature 1 over the rest
# (4 chains, 300 warmup iterations + 300 sampling iterations - use more in practice)
samples_1 <- fit_signatures(counts_21breast, cosmic_signatures_v2[1:4, ],
                            exp_prior = c(10, 1, 1, 1), iter = 600)

# Fit all the signatures, running a single chain for many iterations
# (3000 warmup iterations + 10000 sampling iterations)
samples_2 <- fit_signatures(counts_21breast, cosmic_signatures_v2, chains = 1,
                            iter = 13000, warmup = 3000)
}
}
