# bgchm R package

This is an R package for Bayesian analyses of population genomic data from hybrid zones, including Bayesian genomic cline analysis, estimation of hybrid indexes and ancestry classes, some geographic cline analyses, and accessory plotting functions. This package using Hamiltonian Monte Carlo (HMC) for sampling posterior distributions, with HMC sampling implemented via Stan.

# Installation

You can install this package within R directly from GitHub.

```{R}
## install and load devtools
install.packages("devtools")
library(devtools)
## install bgc-hm
## this will take a bit as it requires compiling a substantial amount of C++ code
devtools::install_github("zgompert/bgc-hm")
## load the bgc-hm package
library(bgchm)
```

# Usage

This is a working version of the software, but not all options have been implemented yet. I am actively developing this package and will post details on usage once I have a version with everything implemented and working (hopefully by the end of Nov. 2023).

# Citations

The general hierarchical Bayesian model used for Bayesidan genomic cline analysis was described here:

The current set of models based on the log-logistic function and using HMC were first described here:

The following paper (in the works) describes this specific R package:



