
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to hemRNA

<!-- badges: start -->
<!-- badges: end -->

`hemRNA` is a free open-source software developed to help clinicians and
biologists in analyzing RNAseq data. It has been designed to
automatically process results in a standardized structure as provided by
the associated script (RNAseq\_main\_script.sh).

`hemRNA` has been developed with the R software, using the [Shiny
package](https://shiny.rstudio.com/).
[Golem](https://github.com/ThinkR-open/golem) has been used for package
compilation and deployment.

## Installation

You can install the development version from
[GitHub](https://github.com/VincentAlcazer/hemRNA) by either cloning the
repository or directly downloading the package in R:

``` r
## Github install 
 install.packages("remotes")
 remotes::install_github("VincentAlcazer/hemRNA", force = T)
 
 hemRNA::run_app()
 
## Local file install
# Download and unzip hemRNA folder from github, then run from R:
 
install.packages("path_to_hemRNA-master_folder", repos = NULL, type="source", force = T)
 
```
