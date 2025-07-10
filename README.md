
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VacBioTLF

<!-- badges: start -->
<!-- badges: end -->

The goal of VacBioTLF is to …

## Installation

You can install the development version of VacBioTLF from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("abby529/VacBioTLF")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(VacBioTLF)
#> Loading required package: dplyr
#> Warning: package 'dplyr' was built under R version 4.3.3
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: tidyr
#> Warning: package 'tidyr' was built under R version 4.3.3
#> Loading required package: stringr
#> Warning: package 'stringr' was built under R version 4.3.3
#> Loading required package: ggpubr
#> Warning: package 'ggpubr' was built under R version 4.3.3
#> Loading required package: ggplot2
#> Loading required package: reporter
#> Warning: package 'reporter' was built under R version 4.3.3
#> Loading required package: common
#> Warning: package 'common' was built under R version 4.3.3
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
