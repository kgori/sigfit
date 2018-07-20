# sigfit
### Flexible Bayesian inference of mutational signatures

sigfit is an R package to estimate signatures of mutational processes and their activities on mutation count data. Starting from a set of single-nucleotide variants (SNVs), it allows both estimation of the exposure of samples to predefined mutational signatures (including whether the signatures are present at all), and identification signatures _de novo_ from the mutation counts. These two procedures are often called, respectively, signature fitting and signature extraction. Furthermore, the signature fitting and extraction methods in sigfit can be seamlessly applied to mutational profiles beyond SNV data, including insertion/deletion (indel) or rearrangement count data. The package provides a range of functions to generate publication-quality graphics of the corresponding mutational catalogues, signatures and exposures.


## Installation
sigfit is an R package. As it is in early development it is not yet on CRAN, but can be installed from inside an R session using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

    devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)


## Usage guide

See the package vignette for detailed usage examples:

    browseVignettes("sigfit")

You can also browse the package vignette in [GitHub](http://htmlpreview.github.io/?https://github.com/kgori/sigfit/blob/master/inst/doc/sigfit_vignette.html).


## Citation

To cite sigfit in publications, please use:

* **Kevin Gori, Adrian Baez-Ortega. sigfit: flexible Bayesian inference of mutational signatures. _bioRxiv_ (2018). doi: 10.1101/372896.**

The corresponding BibTeX entry is:

    @Article{,
        title = {sigfit: flexible Bayesian inference of mutational signatures},
        author = {Kevin Gori and Adrian Baez-Ortega},
        journal = {bioRxiv},
        year = {2018},
        doi = {10.1101/372896},
    }


## Licence

Authors: Kevin Gori and Adrian Baez-Ortega  
Transmissible Cancer Group, University of Cambridge

sigfit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
