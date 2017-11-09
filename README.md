# sigfit
### Discovering mutational signatures through Bayesian inference 

sigfit is an R package to estimate signatures of mutational processes and their activities on mutation count data. Starting from a set of single-nucleotide variants (SNVs), it allows both estimation of the exposure of samples to predefined mutational signatures (including whether the signatures are present at all), and identification signatures _de novo_ from the mutation counts. These two procedures are often called, respectively, signature fitting and signature extraction. Furthermore, the signature fitting and extraction methods in sigfit can be seamlessly applied to mutational profiles beyond SNV data, including insertion/deletion (indel) or rearrangement count data. The package provides a range of functions to generate publication-quality graphics of the corresponding mutational catalogues, signatures and exposures.

## Installation
sigfit is an R package. As it is in early development it is not yet on CRAN, but can be installed from inside an R session using the devtools library.

    devtools::install_github("kgori/sigfit", build_vignettes = TRUE)
    

## Usage guide

See the package vignette for detailed usage examples:

    browseVignettes("sigfit")

