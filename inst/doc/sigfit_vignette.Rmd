---
title: "Fitting and extracting mutational signatures with sigfit"
author: "Kevin Gori and Adrian Baez-Ortega (2018)"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Usage guide}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library(sigfit))
par(mar = c(6, 4, 6, 4))
```

## Introduction

sigfit is used to estimate signatures of mutational processes and their degree of activity on a collection of cancer (or normal) samples. Starting from a set of single-nucleotide variants (SNVs), it allows both estimation of the exposure of samples to predefined mutational signatures (including whether the signatures are present at all), and identification signatures _de novo_ from the mutation counts. These two procedures are often called, respectively, signature fitting and signature extraction. In addition, sigfit includes 'Fit-Ext' models that allow combining signature fitting and extraction in a single inference process. Furthermore, the signature analysis methods in sigfit can be seamlessly applied to mutational profiles beyond SNV data, including insertion/deletion (indel) or rearrangement count data. The package also provides a range of functions to generate publication-quality graphics of the corresponding mutational catalogues, signatures and exposures.

## Installation

sigfit is an R package. As it is in early development it is not yet on CRAN, but can be installed from GitHub using the devtools package.

```{r devtools_instructions, eval=FALSE}
devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)
```

## Usage guide

## Example 1: Fitting mutational signatures to a single simulated sample

This example will use the mutational signatures from [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures) to generate simulated mutation counts, and then use sigfit to fit the signatures back to the simulated data.

First of all we need some mutational signatures to fit to our data. The line below loads the mutational signatures published in COSMIC.

```{r fetch}
data("cosmic_signatures", package = "sigfit")
```

Let's use these signatures to simulate some mutation data. This code will generate 20,000 mutations from a 4:3:2:1 mixture of signatures 1, 3, 7 and 11.

```{r sim}
set.seed(1)
probs <- c(0.4, 0.3, 0.2, 0.1) %*% as.matrix(cosmic_signatures[c(1, 3, 7, 11), ])
mutations <- matrix(rmultinom(1, 20000, probs), nrow = 1)
colnames(mutations) <- colnames(cosmic_signatures)
```

Here is what our simulated counts look like:
```{r plotsim, fig.width=17, fig.height=6.5, out.width="100%", echo=-1}
par(mar = c(4.5,5.5,6.5,1))
sigfit::plot_spectrum(mutations, name = "Simulated counts")
```

### Fitting mutational signatures

Next, we can estimate the exposure of the data to each signature (pretending we ignore that it was generated from signatures 1, 3, 7 and 11). sigfit uses [Stan](http://mc-stan.org/) to run a Bayesian model that produces Markov chain Monte Carlo (MCMC) samples. Arguments to the ```rstan::sampling``` function, such as ```iter```, ```warmup```, etc., can be passed through. For further sampling options, type ```?rstan::sampling``` to read the documentation.

__In general, one should run as many MCMC iterations (```iter``` argument) as one's computer and patience allow, runtime being the major constraint.__

We recommend that the number of warmup (burn-in) iterations (```warmup``` argument) be between one-third and half the value of ```iter```. The behaviour of the MCMC sampler (which ultimately affects the quality of the analysis) depends on parameters set during the warmup, so it is important to run plenty of warmup iterations. By default, ```rstan::sampling``` uses ```iter = 2000``` and ```warmup = floor(iter/2)```; we do not normally recommend going below these values. The ```seed``` argument can be used to make the MCMC samples reproducible over different runs.

We can use ```fit_signatures``` to fit the COSMIC signatures to the simulated counts as follows.

```{r fitting, warning=FALSE}
mcmc_samples_fit <- sigfit::fit_signatures(counts = mutations, 
                                           signatures = cosmic_signatures,
                                           iter = 2000, 
                                           warmup = 1000, 
                                           chains = 1, 
                                           seed = 1)
```

### Retrieving signature exposures

Once we have the result of the MCMC sampling in ```mcmc_samples_fit```, we can retrieve the estimated exposures from it using the ```retrieve_pars``` function. This returns a named list with three matrices, one containing the mean exposures, and the others containing the values corresponding to the lower and upper limits of the highest posterior density (HPD) interval (the Bayesian alternative to a confidence interval) for each exposure in each sample. The ```prob``` argument can be used to indicate the target probability content of the HPD interval (by default, 95% HPD intervals are returned).

Since we are fitting known signatures and not extracting new ones, the exposures will be automatically labelled with the signature names. If the signatures have no names (or have been extracted _de novo_ instead of fitted, as we will see below), they will be labelled as 'Signature A', 'Signature B', etc.

```{r retrieve_exp}
exposures <- sigfit::retrieve_pars(mcmc_samples_fit, 
                                   par = "exposures", 
                                   hpd_prob = 0.90)
names(exposures)
exposures$mean
```

The entire posterior distribution of the signature exposures and other model parameters in the ```mcmc_samples_fit``` object can be further explored by means of the functions provided by the [rstan](http://mc-stan.org/rstan) package. In addition, [ShinyStan](http://mc-stan.org/users/interfaces/shinystan) can be easily used in R for visual exploration of the MCMC samples.

### Visualisation

sigfit provides several easy-to-use plotting functions. As seen in the previous section, the ```plot_spectrum``` function allows visualisation of both mutational catalogues and mutational signatures. These plots can be produced even for catalogues and signatures with arbitrary mutation types (for example, indel or rearrangement signatures).

The ```plot_exposures``` function produces a barplot of the estimated signature exposures in each sample. It needs to be supplied with either the object resulting from MCMC sampling (```mcmc_samples``` argument) or the exposures themselves (```exposures``` argument), the latter being either a matrix, or a list like the one returned by the ```retrieve_pars``` function (above). In the present case, since we have the stanfit object generated by ```fit_signatures```, we will make use of the ```mcmc_samples``` argument.

```{r plot_exp, fig.width=17, fig.height=7, out.width='100%', fig.align="center", echo=-1}
par(mar=c(8,5,3.5,0))
sigfit::plot_exposures(mcmc_samples = mcmc_samples_fit)
```

The bars in this plot are coloured blue if the estimated exposure value is 'sufficiently non-zero'. It is difficult for the model to make hard assignments of which signatures are present or absent due to the non-negative constraint on the estimate, which means that the range of values in the sample will not normally include zero. In practice, 'sufficiently non-zero' means that the lower end of the Bayesian HPD credible interval is above a threshold value close to zero (by default 0.01, and adjustable via the ```thresh``` argument). In this example, sigfit has identified the 4 signatures used to construct the sample.

Next, we would recommend running ```fit_signatures``` again, this time to fit only those signatures (i.e. those rows of the ```cosmic_signatures``` matrix) which have been highlighted as 'sufficiently non-zero' in the plot above, in order to obtain more accurate estimates. We will skip this step in the present example.

We can also examine how effectively the estimated signatures and/or exposures can reconstruct the original count data, using the ```plot_reconstruction``` function. 

__Note that the plotting functions in sigfit are designed with a preference for plotting directly to a PDF file__.

The path to our desired output PDF can be provided using the ```pdf_path``` argument, and each function will automatically select the most appropriate size and graphical parameters for the plot. (We will not make use of this option in the present example, however.) The ```sig_color_palette``` argument can be used to specify custom colours for the signatures in the reconstructed spectrum.

```{r reconstruct, fig.width=25, fig.height=18.5, out.width='100%', warning=FALSE, results="hide", echo=-1}
par(mar=c(5,6,6.5,1))
sigfit::plot_reconstruction(mcmc_samples = mcmc_samples_fit,
                            pdf_path = NULL)
```

The ```plot_all``` function provides a simple way of combining the ```plot_spectrum```, ```plot_exposures``` and ```plot_reconstructions``` functions. The ```plot_all``` function shares most arguments with the other plotting functions, and is useful to avoid running all the plotting functions individually. This function plots exclusively to PDF files, and the ```out_path``` argument is used to indicate the path of the directory where the files should be created (if the directory does not yet exist, it will be automatically created). The ```prefix``` argument applies to the output file names, and can be used to distinguish different sets of plots from each other.

```{r plot_all, eval=FALSE}
sigfit::plot_all(mcmc_samples = mcmc_samples_fit, 
                 out_path = "your/output/dir/here",
                 prefix = "Fitting")
```


## Example 2: Extracting signatures from multiple breast cancer samples

### Generating mutational catalogues

In this second example, we will use single-nucleotide variant (SNV) data from the set of 21 breast cancer samples presented by [Nik-Zainal _et al._ (2012)](http://dx.doi.org/10.1016/j.cell.2012.04.024). These data can be accessed using ```data("variants_21breast")```.

```{r load_mutations}
data("variants_21breast", package = "sigfit")
head(variants_21breast)
```

This table illustrates the structure of the variant data that can be used as input for the package (unless you already have mutational catalogues for your samples). It is a matrix or data frame with one row per variant, and four (or five) columns:

* __Sample ID__ (character, e.g. ```"Sample 1"```).
* __Reference base__ (character: ```"A"```, ```"C"```, ```"G"```, or ```"T"```).
* __Mutated base__ (character: ```"A"```, ```"C"```, ```"G"```, or ```"T"```).
* __Trinucleotide context__ of the variant (character; reference sequence between the positions immediately before (-1) and after (+1) the variant, e.g. ```"TCA"```). This can be obtained from the reference genome that was used to call the variants, using an R package like [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html); however, sequence context information is sometimes provided by variant callers within the 'INFO' field of the VCF file.
* __Optionally__: if information is available about the __transcriptional strand__ where each mutation occurred, this can be incorporated as a fifth column taking character values ```"1"``` / ```"U"``` (for the untranscribed strand) or ```"-1"``` / ```"T"``` (for the transcribed strand). If the table is a data frame, integer values ```1``` and ```-1``` are also admitted. If this column is present in the table, all the estimation and plotting functions will automatically incorporate transcriptional strand information. The interpretation of ```1``` and ```-1``` as the untranscribed and transcribed strands, respectively, allows direct use of the 'Strand' field that is normally provided by variant annotation tools, such as the [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html). (Annotation tools use 'Strand' values of ```1``` and ```-1``` to indicate that the corresponding gene is located on the forward or reverse strand, respectively. Since mutations are always called on the forward strand, this value also indicates the transcriptional strand of mutations.)

Importantly, since a variant can only have a single sample ID, variants which are found in more than one sample need to be included multiple times in the table, using different sample IDs. The order in which the samples are found in this table is the order in which they will be displayed thereafter. In this case, the samples are already sorted alphabetically:
 
```{r show_samples}
unique(variants_21breast[, 1])
```

The first step is to transform these variants into mutational catalogues, which is done by the ```build_catalogues``` function. (You can skip this step if you already have mutational catalogues for each of your samples.)

```{r build_catalogues}
counts_21breast <- sigfit::build_catalogues(variants_21breast)
dim(counts_21breast)
counts_21breast[1:5, 1:8]
```

The mutational catalogues are stored as a matrix of mutation counts, where each row refers to a sample and each column corresponds to a trinucleotide mutation type.

(This example set of 21 mutational catalogues can also be loaded directly using ```data("counts_21breast", package = "sigfit")```).

We can plot the spectrum of all the mutational catalogues using the ```plot_spectrum``` function, as in the previous example. For tables containing more than one catalogue, this function will produce one plot per catalogue, which makes using an output PDF file (```pdf_path``` argument) more convenient. In this example, however, we will plot all the catalogues together.

```{r plot_spectra, fig.width=22, fig.height=23, out.width='100%', fig.align="center", echo=-1}
par(mar = c(5,6,7,2))
par(mfrow = c(7, 3))
sigfit::plot_spectrum(counts_21breast)
```

### Extracting mutational signatures _de novo_

To extract signatures from this set of catalogues, we use the ```extract_signatures``` function, specifying the number of signatures to extract via the ```nsignatures``` argument; this can be a single integer or a range, e.g. ```3:6```. Our recommended approach is first running the function for a small number of iterations and a reasonably wide range of numbers of signatures (e.g. ```nsignatures = 2:8```). When ```nsignatures``` is a range of values, sigfit will automatically determine the most plausible number of signatures present in the data (which is done by assessing goodness of fit through the ```plot_gof``` function). Also, for ranges of ```nsignatures``` the result will be a list, where element ```[[N]]``` in the list corresponds to the extraction results for ```nsignatures = N```.

```{r extraction, eval=FALSE}
mcmc_samples_extr <- sigfit::extract_signatures(counts = counts_21breast,
                                                nsignatures = 2:7,
                                                iter = 1000, 
                                                seed = 1)
```

```{r plot_gof_silent, echo=FALSE, fig.width=9, fig.height=6, out.width="100%"}
## Plot precalculated GOF in order to avoid running the model
data("sigfit_vignette_data", package = "sigfit")
plot(nS, gof, type = "o", lty = 3, pch = 16, col = "dodgerblue4",
     main = paste0("Goodness of fit (", stat, ")\nmodel: NMF"),
     xlab = "Number of signatures", 
     ylab = paste0("Goodness of fit (", stat, ")"))
points(nS[best], gof[best], pch = 16, col = "orangered", cex = 1.1)
cat("Estimated best number of signatures:", nS[best], "\n")
```

The plot above shows that the most plausible number of signatures is four, based on the evolution of the goodness of fit (reconstruction accuracy measured through cosine similarity).

Next, we would recommend running ```extract_signatures``` again, this time with ```nsignatures = 4``` and a much greater number of iterations, in order to obtain more accurate estimates. We will skip this step in the present example.

### Retrieving and plotting extracted signatures

As in the case of signature fitting (Example 1 above), the extracted signatures and exposures can be retrieved using the ```retrieve_pars``` function with ```par = "signatures"``` or ```par = "exposures"```. However, because the output is a list containing the results for each value in ```nsignatures```, __we need to specify which  number of signatures we are interested in, using the ```[[ ]]``` operator__. Note that this is not necessary if a single value of ```nsignatures``` was used for extraction.

```{r extr_names, eval=FALSE}
names(mcmc_samples_extr)
```
```{r extr_names_silent, echo=FALSE}
print(c("nsignatures=1", "nsignatures=2", "nsignatures=3", "nsignatures=4", "nsignatures=5", "nsignatures=6", "nsignatures=7", "best"))
```
```{r retrieve_sigs, eval=FALSE}
## Note: mcmc_samples_extr[[N]] contains the extraction results for N signatures
extr_signatures <- sigfit::retrieve_pars(mcmc_samples_extr[[4]],
                                         par = "signatures")
```
```{r show_signames}
rownames(extr_signatures$mean)
```

Plotting can be done through the functions seen in Example 1, with the difference that, when plotting directly from MCMC results (via the ```mcmc_samples``` argument), we need to select the relevant element of the results list using ```mcmc_samples_extr[[N]]```, as seen above. Note that this is not necessary if a single value of ```nsignatures``` was used for extraction.

Below we plot the signatures extracted from these 21 catalogues.

```{r plot_sigs, warning=FALSE, fig.width=25, fig.height=10, out.width='100%', fig.align="center", echo=-1}
par(mar = c(6,7,6,1))
par(mfrow = c(2, 2))
sigfit::plot_spectrum(extr_signatures)
```

These are a combination of COSMIC signatures 1, 2, 3, 5 and 13. Note that the signatures published in [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures) were obtained using a collection of hundreds of catalogues across many cancer types, which offered much higher statistical power than the 21 breast cancer catalogues employed here. Furthermore, the signatures obtained by sigfit show high similarity to those originally reported by [Nik-Zainal _et al._ (2012)](http://dx.doi.org/10.1016/j.cell.2012.04.024) (Fig. 2A). Note that signatures C and D in Nik-Zainal _et al._, which are very similar, have been identified by sigfit as a single signature (Signature B in the plot above).

### Using the EMu (Poisson) signature model

By default, both ```fit_signatures``` and ```extract_signatures``` make use of an NMF-inspired Dirichlet–Multinomial model of signatures, which is statistically equivalent to the non-negative matrix factorisation (NMF) approach pioneered by [Alexandrov _et al._ (2013)](https://www.nature.com/articles/nature12477). Alternatively, users who are interested in the Poisson model presented by [Fischer _et al._ (2013)](https://doi.org/10.1186/gb-2013-14-4-r39) can use the ```model = "emu"``` option to select this model, which is able to account for variation in mutational opportunity (the opportunity for each mutation type to occur in each sample's genome; this is specified via the ```opportunities``` argument). For further details, type ```?extract_signatures``` to read the documentation.

Although signature representations differ between the NMF model and the EMu model (insofar as signatures obtained through the latter are not relative to the mutational opportunities of a specific genome/exome), signatures can be converted between both model representations by means of the ```convert_signatures``` function. For further details, type ```?convert_signatures``` to read the documentation.

### Converting between genome- and exome-relative representations of signatures

In some analyses, signatures that have been inferred from whole-genome mutation data need to be fitted to mutation data from whole-exome sequencing samples. For instance, if the cohort of samples includes both whole genomes and exomes, it is normally advisable to extract signatures from the whole-genome samples (provided that there are enough of them) and re-fit the resulting signatures to the whole-exome samples. (An alternative approach would be using the EMu/Poisson model mentioned above, with a matrix of sample-specific mutational opportunities that reflects the origin of each catalogue.)

To convert a matrix of genome-derived mutational signatures into their exome-relative representations, the ```convert_signatures``` function needs to be applied twice: once to normalise the signatures by the whole-genome mutational opportunities (using arguments ```ref_opportunities = "human-genome", model_to = "emu"```), and once again to impose whole-exome opportunities (using ```ref_opportunities = "human-exome", model_to = "nmf"```). An example is shown below.

```{r convert_sigs, eval=F}
# Following from the signature extraction example above
genome_signatures <- sigfit::retrieve_pars(mcmc_samples_extr[[4]],
                                           par = "signatures")
normalised_signatures <- sigfit::convert_signatures(genome_signatures, 
                                                    ref_opportunities = "human-genome",
                                                    model_to = "emu")
exome_signatures <- sigfit::convert_signatures(normalised_signatures, 
                                               ref_opportunities = "human-exome",
                                               model_to = "nmf")

par(mfrow = c(2, 1))
sigfit::plot_spectrum(genome_signatures$mean[4,], name = "Signature D, Genome-relative probabilities")
sigfit::plot_spectrum(exome_signatures[4,], name = "Signature D, Exome-relative probabilities")
```

```{r convert_sigs_silent, echo=F, fig.width=20, fig.height=14, out.width="100%"}
genome_signatures <- extr_signatures
normalised_signatures <- sigfit::convert_signatures(genome_signatures, ref_opportunities="human-genome", model_to="emu")
exome_signatures <- sigfit::convert_signatures(normalised_signatures, ref_opportunities="human-exome", model_to="nmf")
par(mfrow = c(2, 1), mar = c(5,5.5,6.5,1))
sigfit::plot_spectrum(genome_signatures$mean[4,], name="Signature D, Genome-relative probabilities")
sigfit::plot_spectrum(exome_signatures[4,], name="Signature D, Exome-relative probabilities")
```

### Using Fit-Ext models to discover rare or novel signatures

One novelty in sigfit is the introduction of Fit-Ext models, which are able to extract mutational signatures _de novo_ while fitting a set of predefined signatures that are already known to be present in the data. Such models are useful for the discovery of rare signatures for which there is some qualitative evidence, but insufficient support as to allow deconvolution via conventional methods; or in cases where only signature fitting is possible, yet the data clearly display a mutational pattern which cannot be captured by the available signatures. For further details on the Fit-Ext models, please refer to the [sigfit paper](https://www.biorxiv.org/content/early/2018/07/20/372896).

The Fit-Ext models can be accessed via the ```fit_extract_signatures``` function. This is used similarly to ```extract_signatures```, with the exception that a matrix of known signatures to be fitted needs to be provided via the ```signatures``` argument (as in ```fit_signatures```), and that the number of additional signatures to extract is provided via the ```num_extra_sigs``` argument. Unlike the ```nsignatures``` argument in ```extract_signatures```, ```num_extra_sigs``` currently admits only scalar values. For further details, type ```?fit_extract_signatures``` to read the documentation.

___

sigfit is an R package developed by the [Transmissible Cancer Group](http://www.tcg.vet.cam.ac.uk/) in the University of Cambridge Department of Veterinary Medicine.
