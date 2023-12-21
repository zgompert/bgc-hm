# bgchm R package

This is an R package for Bayesian analyses of population genomic data from hybrid zones, including Bayesian genomic cline analysis, estimation of hybrid indexes and ancestry classes, some geographic cline analyses, and accessory plotting functions. This package using Hamiltonian Monte Carlo (HMC) for sampling posterior distributions, with HMC sampling implemented via Stan.

# Installation

You can install this package within R directly from GitHub.

```R
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

# Examples

Fit genomic clines for an example data set with known genotypes. This data set comprises 51 ancestry-informative loci, 100 putative hybrids and 50 individuals represntative of each parental species. All loci are diploid. The data were simulated with dfuse using an underdominance model with xxx (the underdominance model is described in [Fierno et al. 2023](https://cshperspectives.cshlp.org/content/early/2023/09/21/cshperspect.a041434).
```R
## load the data set
data(genotypes)
## this includes three objects, GenHybrids, GenP0, and GenP1

## estimate parental allele frequencies, uses default HMC settings
p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")

## estimate hybrid indexes, uses default HMC settings
## and uses point estimates (posterior medians) of allele frequencies
h_out<-est_hi(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")

## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## this too uses default HMC settings
gc_out<-est_genocl(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="genotype",ploidy="diploid",hier=TRUE)

```

# Citations

The general hierarchical Bayesian model used for Bayesidan genomic cline analysis was described here:

[Gompert Z, Buerkle CA (2011) Bayesian estimation of genomic clines. Molecular Ecology, 20:2111-2127.](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2011.05074.x)

The current set of models based on the log-logistic function and using HMC were first described here:

[Fierno TJ, Semenov G, Dopman EB, Taylor SA, Larson EL, Gompert Z (2023) Quantitative analyses of coupling in hybrid zones. Cold Spring Harb Perspect Biol, a041434.](https://cshperspectives.cshlp.org/content/early/2023/09/21/cshperspect.a041434)

The following paper (in the works) describes this specific R package:



