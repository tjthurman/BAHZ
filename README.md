
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![Lfecycle experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) ![CRAN status](http://www.r-pkg.org/badges/version/bahz) <!-- badges: end -->

bahz
====

**bahz** (**B**ayesian **a**nalysis of **h**ybrid **z**ones) is an R package for fitting and analyzing one-dimensional, geographic cline models in hybrd zones.

**bahz** is currently under active development.

Installation
------------

Before installing **bahz**, you'll need to install the R package RStan, which interfaces with the Stan modelling language. To install RStan, follow the instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

Once Rstan is installed, you can install **bahz** from Github. If you don't already have it, you'll need to install the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
install.packages("devtools")
```

Then, you can install **bahz**:

``` r
library(devtools)
install_github("tjthurman/BAHZ", ref = "master") 
library(blirp)
```

**bahz** is not currently available from CRAN. If you have any issues installing **bahz**, please [file an issue on GitHub](https://github.com/tjthurman/BAHZ/issues).

Citation
--------

No citation yet.

License
-------

**bahz** is distributed for free under the [GNU Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html).
