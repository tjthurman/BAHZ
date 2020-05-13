# bahz 0.0.0.9018

New features/changes:

* Restructured package according to updated recommended structure in `rstantools` version `2.0.0`.


# bahz 0.0.0.9017

New features/changes:

* `cline_summary()` now provides both posterior means and posterior medians for each parameter.

* `predict_cline()` outputs best fit estimates based on both posterior means and posterior medians, and `plot_cline()` can plot best-fit lines based on either of these best-fit estimates (defaults to posterior means).

* `predict_cline()` and `plot_cline()` can now output both HPDI credible intervals (the default) and equal-tail (ET) credible intervals.

# bahz 0.0.0.9016

Changes:

* Deleting unnecessary files that prevented installation on Windows and Linux.

* Remake corrupted file in test suite.

# bahz 0.0.0.9015

Changes:

* Consolidated `plot_geno_cline()` and `plot_pheno_cline()` into a single function, `plot_cline()`.

Behind the scenes:

* Streamlining of testing and argument checking for `plot_cline()`, `predict_cline()`, and `cline_summary()` functions. 

# bahz 0.0.0.9014

New features:

* Ability to add credible intervals when predicting clines and plotting clines.

* Updated phenotypic models to allow three different variance structures: constant across sites, differing across sites, or a pooled hierachical model. 


# bahz 0.0.0.9013

Bugfix update

* Fixes #57- A bug which caused bahz to always estimate increasing phenotypic clines, even when data were for a decreasing cline.

# bahz 0.0.0.9012

Major update implementing phenotypic clines:.

* Fully implementeded for phenotypic clines (of normally distributed phenotypes), from simulated data to plotting. 

* Improved template for prior specifications with more explanation.

* Unifies predict_*_cline functions.

* Adds support for new prior distributions: beta distributions for pmin, pmax, tauL, rauR, tauM, and f, and normal distributions for pmin and pmax. 

* Added a citation file.

* Some minor typo fixes. 

# bahz 0.0.0.9011

First entry in the NEWS.md file. Beta, working version of the package. 

* Fully implemented for genetic clines: from simulated data to plotting

* All functions documented

* All functions tested

* R CMD check passes with 1 note (on GNU make requirement). 

