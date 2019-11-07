<p align="center"><img src="logo.png" alt="sigfit" width="700"/></p>

## Flexible Bayesian inference of mutational signatures

[![Build Status](https://travis-ci.org/kgori/sigfit.svg?branch=master)](https://travis-ci.org/kgori/sigfit)


__sigfit__ is an R package to estimate signatures of mutational processes and their activities on mutation count data. Starting from a set of single-nucleotide variants (SNVs), it allows both estimation of the exposure of samples to predefined mutational signatures (including whether the signatures are present at all), and identification of signatures _de novo_ from the mutation counts. These two procedures are often called, respectively, signature fitting and signature extraction. In addition, sigfit implements novel methods to combine signature fitting and extraction in a single inferential process, thus facilitating the deconvolution of rare or admixed signatures. The package provides interfaces to four different Bayesian models of signatures (multinomial, Poisson, normal and negative binomial), as well as a range of functions to generate publication-quality graphics of the corresponding mutational catalogues, signatures and exposures. Furthermore, the signature fitting and extraction methods in sigfit can be seamlessly applied to mutational profiles beyond SNV data, including indel or rearrangement count data, and even real-valued data such as DNA methylation profiles.

### Enhancements in version 2.0 (November 2019)

* New models for analysis of real-valued data ("normal") and robust fitting to sparse data ("negbin")
* New COSMIC v3 signatures (67 SBS signatures) and test datasets
* Straightforward analysis and plotting of signatures defined over arbitrary mutation types
* Support for mutational opportunities in all signature models
* Support for signature and exposure priors in all signature models
* Enhanced plotting functionalities
* Increased MCMC sampling efficiency
* Extended package vignette

## Installation
sigfit is an R package. As it is in early development it is not yet on CRAN, but can be installed from inside an R session using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

    devtools::install_github("kgori/sigfit", build_vignettes = TRUE,
                             build_opts = c("--no-resave-data", "--no-manual"))
    
The arguments `build_vignettes` and `build_opts` are necessary for the package vignette to be built.

For solutions to some of the problems that may arise during installation, see the __[Troubleshooting installation](#troubleshooting-installation)__ section.


## Usage guide

See the package vignette for detailed usage examples:

    browseVignettes("sigfit")
    
If the vignette is not available because it was not built during installation, you can download it and build it using the `rmarkdown` package (the vignette will be saved as an HTML file in your current working directory):

    download.file("https://raw.githubusercontent.com/kgori/sigfit/master/vignettes/sigfit_vignette.Rmd",
                  destfile = "sigfit_vignette.Rmd")
    rmarkdown::render("sigfit_vignette.Rmd")
    browseURL("sigfit_vignette.html")

### You can also browse the package vignette on [GitHub](http://htmlpreview.github.io/?https://github.com/kgori/sigfit/blob/master/doc/sigfit_vignette.html).


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


## Troubleshooting installation

Below are the solutions to some of the problems that may arise during installation.

__Problem:__

    Error: 'rstan_config' is not an exported object from 'namespace:rstantools'
    
__Solution:__  
Update rstantools: `devtools::install_github("stan-dev/rstantools")`

---

__Problem:__

    C++14 standard requested but CXX14 is not defined
    
__Solution:__  
Provide R with c++14 options via the file `~/.R/Makevars`, e.g.

    CXX14 = g++
    CXX14FLAGS = -g -O2
    CXX14PICFLAGS = -fpic
    CXX14STD = -std=gnu++14

---

__Problem:__

    make: *** [stanExports_sigfit_ext.o] Error 1
    ERROR: compilation failed for package ‘sigfit’

Which is preceded by a series of compilation errors related to `StanHeaders`, such as:

    .../StanHeaders/include/stan/math/rev/mat/functor/adj_jac_apply.hpp:619:15: error: invalid use of ‘auto’

__Solution:__  
Provide R with c++14 options via the file `~/.R/Makevars`, as described above.

---

__Problem:__

    g++: error: unrecognized command line option '-std=gnu++14'

__Solution:__
Upgrade the `gcc` and `g++` compilers to version 5. In __Ubuntu__, this can be done as follows:

    sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    sudo apt-get update
    sudo apt-get install gcc-5 g++-5
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5

Alternatively, if you have version __4.8.1__ or higher of `gcc` and `g++`, you can enable C++14 features by substituting `-std=gnu++1y` for `-std=gnu++14` in `~/.R/Makevars` (see above).
