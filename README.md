
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tradisanalyseR

<!-- badges: start -->
<!-- badges: end -->

tradisanalyseR was developed as a collection of scripts to further
analyse output produced by the Bio-TraDIS toolkit
(<https://github.com/sanger-pathogens/Bio-Tradis>). It takes the
resulting csv files for any number of conditions and formats it for use
in these scripts. This is a work-in-progress. Please email Geri at
<geraldine.sullivan@hdr.mq.edu.au> if you have any reproducibility
issues.

## Scripts
##### structure_csv()
Takes all csv files from a folder and combines them into one data table for visualisation and further processing.

##### logfc_plots()
Using the structure_csv() output and a gene list, plots log fold changes for those genes over all conditions in one graph.

##### diagnostics_rc()
Takes all sites.csv files in a folder and plots a PCA to determine whether batch effects are present.

##### diagnostics_fc()
Takes all csv output files from tradis_comparison.R in a folder and plots fold changes along the chromosome to determine whether chromosomal bias is present. Currently plots p-values as well but will be removed.

##### get_operons()
Requires EcoCyc/BioCyc output of transcription units, and plots all significant genes and their operons to determine trends.

##### get_sig()
Takes all csv output files and produces new output showing only significant genes.

##### promoter_weblogo()
Requires upstream 50bp from all genes (available with Artemis). Takes all significant genes in a condition and plots a sequence logo to identify possible binding motifs for transcription factors.

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
EDASeq (>= 2.24.0)
ggplot2
ggseqlogo
gridExtra
msa
purrr
reshape2
RUVSeq
RWebLogo
seqinr
stats
stringr
```
