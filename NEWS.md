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

