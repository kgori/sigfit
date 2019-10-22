<p align="center"><img src="logo.png" alt="sigfit" width="700"/></p>

## Flexible Bayesian inference of mutational signatures

[![Build Status](https://travis-ci.org/kgori/sigfit.svg?branch=master)](https://travis-ci.org/kgori/sigfit)


sigfit is an R package to estimate signatures of mutational processes and their activities on mutation count data. Starting from a set of single-nucleotide variants (SNVs), it allows both estimation of the exposure of samples to predefined mutational signatures (including whether the signatures are present at all), and identification of signatures _de novo_ from the mutation counts. These two procedures are often called, respectively, signature fitting and signature extraction. In addition, sigfit implements novel methodos to combine signature fitting and extraction into a single inferential process. The package provides interfaces to four different Bayesian models of signatures (multinomial, Poisson, normal and negative binomial), as well as a range of functions to generate publication-quality graphics of the corresponding mutational catalogues, signatures and exposures. Furthermore, the signature fitting and extraction methods in sigfit can be seamlessly applied to mutational profiles beyond SNV data, including indel or rearrangement count data, and even real-valued data such as DNA methylation profiles.

__Improvements in version 2.0__

* New signature models for analysis of real-valued catalogues ("normal") and noise-robust signature fitting ("negbin").
* Support for COSMIC signatures v3 (67 SBS signatures).
* Support for non-standard catalogues and signatures with arbitrary mutation categories, including generic spectrum plots.
* Support for mutational opportunities in all signature models.
* Support for signature and exposure priors.
* Enhanced plotting functionalities.
* Increased sampling efficiency.
* Additional test datasets.

## Installation
sigfit is an R package. As it is in early development it is not yet on CRAN, but can be installed from inside an R session using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

    devtools::install_github("kgori/sigfit", build_opts = c("--no-resave-data", "--no-manual"))
    
`build_opts` can be omitted, but it will not build the package vignettes in this case.

#### Troubleshooting installation

Problem:

    Error: 'rstan_config' is not an exported object from 'namespace:rstantools'
    
Solution:
Update rstantools: `devtools::install_github("stan-dev/rstantools")`

Problem:

    C++14 standard requested but CXX14 is not defined
    
Solution:
Provide R with c++14 options via the file `~/.R/Makevars`, e.g.

    CXX14 = g++
    CXX14FLAGS = -g -O2
    CXX14PICFLAGS = -fpic
    CXX14STD = -std=gnu++14


## Usage guide

See the package vignette for detailed usage examples:

    browseVignettes("sigfit")

**You can also browse the package vignette on [GitHub](http://htmlpreview.github.io/?https://github.com/kgori/sigfit/blob/master/doc/sigfit_vignette.html).**


## Citation

To cite sigfit in publications, please use:

* **Kevin Gori, Adrian Baez-Ortega. sigfit: flexible Bayesian inference of mutational signatures. _bioRxiv_, 372896 (2018). doi: [10.1101/372896](http://doi.org/10.1101/372896).**

The corresponding BibTeX entry is:

    @Article{sigfit,
        title = {sigfit: flexible Bayesian inference of mutational signatures},
        author = {Gori, Kevin and Baez-Ortega, Adrian},
        journal = {bioRxiv},
        year = {2018},
        pages = {372896},
        doi = {10.1101/372896}
    }


## Licence

Authors: Kevin Gori and Adrian Baez-Ortega  
Transmissible Cancer Group, University of Cambridge

sigfit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
