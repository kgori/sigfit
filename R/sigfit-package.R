#' sigfit: flexible Bayesian inference of mutational signatures
#'
#' @description This package implements Bayesian models for fitting and extracting mixtures of
#' mutational signatures from mutation count data. It provides interfaces to four different
#' Bayesian signature models (multinomial, Poisson, normal and negative binomial), as well as
#' auxiliary functions for analysing and plotting resulting data. Markov chain Monte Carlo
#' sampling is performed using Stan.
#'
#' @docType package
#' @name sigfit-package
#' @aliases sigfit
#' @useDynLib sigfit, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Gori K, Baez-Ortega A (2018). sigfit: flexible Bayesian inference of mutational signatures. bioRxiv, 372896. DOI: 10.1101/372896.
#'
NULL
