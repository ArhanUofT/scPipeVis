
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scPipeVis

<!-- badges: start -->

<!-- badges: end -->

The goal of scPipeVis is to …

## Description

Single cell RNA sequencing (scRNA-seq) is a genomic approach used to
quantify messenger RNA sequences (whole transcriptome) at the resolution
of a single cell, which can be scaled to tens of thousands of cells at
the same time. This results in generation of a cell by gene counts
matrix (in R gene by cell) as the primary data, and additional metadata
from any given scRNA-seq experiment. Data transformation, modelling, and
analysis methods operate downstream on the above described data.
`scPipeVis` is an R package that lets you perform all the steps with a
single function so users can download the package and run this function
to do a very preliminary bare-bones data analysis of scRNA-seq data.
Additionally there is a visualisation function that creates basic plots
which are typically used for the visualisation of scRNA-seq data. The
complete motivation behind the creation of this package follows. There
are major limitations of this package which are also talked about below.

### Motivation

A default single cell RNA sequencing data analysis pipeline comprises of
several steps, from quality control, normalization of count matrix,
feature selection, dimensionality reduction, clustering, annotating cell
types, and then differential gene expression followed by gene set
enrichment analysis. In my literature review, I found that all of these
steps can be performed in so many different ways by using different
methods at each step. This process is further complicated by tuning of
hyperparameters at many different steps. Therefore, any user while
analyzing their data, has to make so many choices/decisions at each
step. These decisions, like tuning the hyperparameter or choosing the
right method is absolutely not trivial. In other words, these choices at
each step can affect your downstream analyses. Therefore, making well
informed decisions are crucial. This is where scPipeVis comes in,
following an extensive literature review, I have put together a pipeline
that takes care of most of these above mentioned steps using the current
best practices methods based on published reviews benchmarking different
methods for all these steps. Further scPipeVis provides a function to
visualize each step of the pipeline with many different and appropriate
depictions. These depictions are chosen after going through the
literature and figuring out what kind of visualizations are required at
each step. Therefore, this part of the package will facilitate informed
decision making for users as they go on to select the best suitable
methods and hyper parameters which in turn best suits their own datasets
and meets their needs.

### Limitations of `scPipeVis`

This package is created for preliminary data analysis of scRNA-seq data,
at the end of the day the user will have to eventually figure out if
certain methods work better to answer the question that they are trying
to answer using their data. The user will also eventually have to tune
hyperparameters of certain functions, if they dont get satisfactory
results from this pipeline. Last but not the least, here I don’t provide
an exhaustive set of visualizations. I will continue working on the
package to add more and more visualizations and make the pipeline more
end to end.

`scPipeVis` package was developed using `R version 4.5.1 (2025-06-13)`,
`Platform: aarch64-apple-darwin20` and
`Running under: macOS Sonoma 14.3`

## Installation

To install the latest version of scPipeVis:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("ArhanUofT/scPipeVis", build_vignettes = TRUE)
library("scPipeVis")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scPipeVis)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

## Overview

``` r
ls("package:scPipeVis")
#> [1] "run_pipeline"       "visualise_pipeline"
# Fix this: only show your packages
browseVignettes("scPipeVis")
#> No vignettes found by browseVignettes("scPipeVis")
```

`scPipeVis` contains 2 main functions as of now

1.  ***run_pipeline*** for running a typical scRNA-seq pipeline given
    the counts matrix/data in a `SingleCellExperiment` object.

2.  ***visualise_pipeline*** for creating typical plots which are used
    during exploratory scRNA-seq data analysis.

Refer to package vignettes for more details. An overview of the package
is illustrated below.

## Contributions

- `scPipeVis` was developed by Arhan Rupani, 4th Bioinformatics and
  Computational Biology student at the University of Toronto. The author
  did an extensive literature review to understand the field of
  scRNA-seq data analysis and created this pipeline.
- The `run_pipeline` function was written by the author and runs a
  typical scRNA-seq pipeline given the counts matrix/data in a
  `SingleCellExperiment` object. The `visualise_pipeline` fucntion was
  written by the author and creats typical plots which are used during
  exploratory scRNA-seq data analysis. Both the packages make use of
  functions implemented in already existing R packages made for
  scRNA-seq data analysis. These packages are mentioned below and their
  full references are in the References section below.
- TODO Other packages used (just add names)
- Generative AI was not used at any point throughout the process of this
  package’s development

## References

TODO

## Acknowledgements

This package was developed as part of an assessment for 2025 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `scPipeVis` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.
