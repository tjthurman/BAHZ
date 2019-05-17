
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![Lfecycle experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) ![CRAN status](http://www.r-pkg.org/badges/version/bahz) <!-- badges: end -->

bahz
====

bahz (Bayesian Analysis of Hybrid Zones) is an R package for fitting and analyzing one-dimensional, geographic cline models in hybrd zones.

bahz is currently under active development.

Installation
------------

Before installing bahz, you'll need to install the R package RStan, which interfaces with the Stan modelling language. To install RStan, follow the instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

Once Rstan is installed, you can install bahz from Github. If you don't already have it, you'll need to install the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
install.packages("devtools")
```

Then, you can install bahz:

``` r
devtools::install_github("tjthurman/BAHZ") 
library(bahz)
```

bahz is not currently available from CRAN. If you have any issues installing bahz, please [file an issue on GitHub](https://github.com/tjthurman/BAHZ/issues).

How to use bahz
---------------

Though bahz is still under development, the core model-fitting functionality is operational. Here's a mini-tutorial (with minimal explanation) on how to fit a hybrid zone model with bahz. There is further information in the help files for each function. A more complete tutorial is in the works.

#### Prepare Data

To generate data to analyze, we can simulate genetic data from a one-dimenional genetic cline. bahz provides a function for performing this simulation, which automatically outputs the data in the format needed for downstream analyses. See the the help file at `?prep_geno_data` for information on supported formats.

Here, we will simulate a cline that was sampled across a 500 kilometer transect, with sampling sites every 15 kilometers and 30 individuals sampled at each site. For the cline, we will simulate a decreasing cline, with a cline center at 237 km, a width of 62 km, and minimim and maximum allele frequencies of 0.05 and 0.99, respectively:

``` r
set.seed(123)
data <- sim_geno_cline(transect_distances = seq(0,500,15), n_ind = 30, Fis = 0,
                       decrease = T, center = 237, width = 62, 
                       pmin = 0.05, pmax = .99)
head(data)
#>   site transectDist cline.p cline.f AA Aa aa  N emp.p emp.f
#> 1    1            0    0.99       0 30  0  0 30 1.000     0
#> 2    2           15    0.99       0 29  1  0 30 0.983     0
#> 3    3           30    0.99       0 28  2  0 30 0.967     0
#> 4    4           45    0.99       0 30  0  0 30 1.000     0
#> 5    5           60    0.99       0 30  0  0 30 1.000     0
#> 6    6           75    0.99       0 28  2  0 30 0.967     0
```

#### Specify priors

bahz fits Bayesian models, and thus users need to specify prior distributions for the parameters to be estimated. This is done through a prior configuration file. Users first generate a configuration file with the command:

``` r
make_prior_config()
```

By default this will generate a configuration file template, `prior_config_template.yaml`, in the current working directory. Users can specify other file locations or names. Open the file and edit the priors as you see fit. For this tutorial, you can leave the priors at their default vaules.

#### Fit the model

With our data and priors ready, we can now fit a cline model. For now, we will fit the simplest cline model, using the binomial model of allele frequencies and no introgression tails:

``` r
set.seed(456)
cline.fit <- fit_geno_cline(data = data, prior_file = "prior_config_template.yaml",
                            type = "bi", tails = "none")
#> Warning: There were 2080 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: There were 1001 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: There were 2 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> http://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

#### Analyze the results

We can summarize the posterior distribution of the cline parameters:

``` r
cline_summary(cline.fit)
#>    param    mean se_mean      sd low_0.95_HPDI up_0.95_HPDI n_eff Rhat
#> 1 center   15.04   94.54  135.65       -104.41       239.14     2 6.38
#> 2  width   81.15   52.10  138.89          4.21       401.22     7 1.49
#> 3   pmin    0.16    0.04    0.07          0.03         0.20     2 3.26
#> 4   pmax    0.90    0.04    0.06          0.80         0.99     2 4.08
#> 5    dev 2142.53  790.45 1145.03        112.82      3000.82     2 5.43
```

The model has done a good job. The point estimates of the parameters are quite close to the values we specified, and the 95% credible intervals contain the simulated values.

bahz provides some simple plotting functions for visualizating fitted clines:

``` r
plot_geno_cline(stanfit = cline.fit, data = data, add.obs.freqs = T, 
                xlab = "distance", ylab = "allele frequency")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

For users who wish to make higher-quality, customized plots, bahz has a helper function to calculate predicted allele frequencies from fitted clines. Those data can then be used in the plotting system of the user's choice, e.g.:

``` r
pred_cline <- predict_geno_cline(stanfit = cline.fit, data = data)
head(pred_cline)
#> # A tibble: 6 x 2
#>   transectDist     p
#>          <dbl> <dbl>
#> 1        -9    0.727
#> 2        -8.48 0.723
#> 3        -7.96 0.720
#> 4        -7.44 0.716
#> 5        -6.93 0.713
#> 6        -6.41 0.709
```

``` r
# Use the ggplot2 package
library(ggplot2)
ggplot(data = pred_cline, aes(x = transectDist, y = p)) +
  geom_line(color = "orange") +
  xlab("distance along transect") +
  ylab("allele frequency")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Citation
--------

No citation yet.

License
-------

bahz is distributed for free under the [GNU Public License v3](https://www.gnu.org/licenses/gpl-3.0.en.html).
