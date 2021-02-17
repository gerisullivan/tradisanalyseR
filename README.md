
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tradisanalyseR

<!-- badges: start -->

<!-- badges: end -->

tradisanalyseR was developed as a collection of scripts to further
analyse output produced by the Bio-TraDIS toolkit
(<https://github.com/sanger-pathogens/Bio-Tradis>). It takes the
resulting csv files for any number of conditions and the corresponding
EMBL file and formats it for use in these scripts.

## Installation

<!-- You can install the released version of tradisanalyseR from  [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("tradisanalyseR")
```
-->

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gerisullivan/tradisanalyseR")
```

tradisanalyseR has the following dependencies:

``` r
dplyr
ggplot2
purrr
stats
stringr
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tradisanalyseR)
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
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
