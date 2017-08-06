# sigfit
Discovering mutational signatures through Bayesian inference

## Installation
sigfit is an R package. It can be installed from inside an R session using the devtools library

    devtools::install_github("kgori/sigfit", build_vignettes = TRUE)
    

## Usage guide
First of all we need some mutational signatures to fit to our data. This command will fetch mutational signatures from [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures)
    
    sigs <- sigfit::fetch_cosmic_data()

Let's use these signatures to simulate some mutation data.
This code will generate 10000 mutations from a 4:3:2:1 mixture of signatures 1, 5, 7 and 11.

    probs <- as.matrix(sigs[, c(1, 5, 7, 11)]) %*% c(0.4, 0.3, 0.2, 0.1)
    mutations <- as.vector(rmultinom(1, 10000, probs))
    names(mutations) <- rownames(sigs)
    
Now we can estimate the exposure of the data to each signature (pretending we don't already know that
it was generated from 1, 5, 7, 11). ```sigfit``` uses [Stan](http://mc-stan.org/) to run a Bayesian model
that produces Markov Chain Monte Carlo samples.

    mcmc_samples <- sigfit::run_sampling(counts = mutations, signatures = sigs,
                                         iter = 1500, warmup = 500)
                                         
See the vignettes for more a detailed example.
