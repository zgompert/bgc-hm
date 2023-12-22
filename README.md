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

This software now works, but I am still in the process of writing the examples and testing everything. Thus, for the moment use this with caution (and maybe shoot me an email if you use this before I remove this message).

# Examples

**Fit genomic clines for an example data set with known genotypes**. This data set comprises 51 ancestry-informative loci, 100 putative hybrids and 50 individuals represntative of each parental species. All loci are diploid. The data were simulated with dfuse assuming 110 demes, m = 0.1 between neighboring demes, and 10 impacting hybrid fitness via underdominance (the underdominance model for dfuse is described in [Fierno et al. 2023](https://cshperspectives.cshlp.org/content/early/2023/09/21/cshperspect.a041434).
```R
## load the data set
data(genotypes)
## this includes three objects, GenHybrids, GenP0, and GenP1

## estimate parental allele frequencies, uses default HMC settings
p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")

## estimate hybrid indexes, uses default HMC settings
## and uses point estimates (posterior medians) of allele frequencies
h_out<-est_hi(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")

## plot hybrid index estimates with 90% equal-tail probability intervals
## sorted by hybid index, just a nice way to visualize what we have
## in this example few hybrids have intermediate hybrid indexes 
plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:100,h_out$hi[order(h_out$hi[,1]),3],1:100,h_out$hi[order(h_out$hi[,1]),4])

## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
gc_out<-est_genocl(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,n_iters=4000)

## how variable is introgression among loci, lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but think this is mostly a good idea
sz_out<-sum2zero(hmc=gc_out$gencline_hmc,transform=TRUE,ci=0.90)

## plot genomic clines for the 51 loci, first without the sum-to-zero constraint
## then with it... these differ more for some data sets than others
gencline_plot(center=gc_out$center[,1],v=gc_out$gradient,pdf=FALSE)
gencline_plot(center=sz_out$center[,1],v=sz_out$gradient,pdf=FALSE)

## summarize loci with credible deviations from genome-averge gradients, here the focus is
## specifically on steep clines indicative of loci introgressing less than the average
which(sz_out$v[,2] > 1) ## index for loci with credibly steep clines
sum(sz_out$v[,2] > 1) ## number of loci with credibly steep clines

## last, lets look at interspecific ancestry for the same data set, this can
## be especially informative about the types of hybrids present
q_out<-est_Q(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")

## plot the results
tri_plot(hi=q_out$hi[,1],Q10=q_out$Q10[,1],pdf=FALSE,pch=19)
## note that some individuals appear to be likely backcrosses (close to the outer lines of the triangles)
## but the inidividals with intermediate hybrid indexes are clearly not F1s but rather late generation hybrids
```

**Fit genomic clines for an example data set with genotype likelihoods**. This data set comprises 51 ancestry-informative loci, 100 putative hybrids and 50 individuals represntative of each parental species. All loci are diploid. The data were simulated with dfuse assuming 110 demes, m = 0.1 between neighboring demes, and 100 loci weakly impacting hybrid fitness via underdominance.
```R
## load the data set
data(gliks)
## this includes three objects, GlikHybrids, GlikP0, and GlikP1

## estimate parental allele frequencies, uses default HMC settings
p_out<-est_p(G0=GlikP0,G1=GlikP1,model="glik",ploidy="diploid")

## estimate hybrid indexes, uses default HMC settings
## and uses point estimates (posterior medians) of allele frequencies
h_out<-est_hi(Gx=GlikHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="glik",ploidy="diploid")

## plot hybrid index estimates with 90% equal-tail probability intervals
## sorted by hybid index, just a nice way to visualize what we have
## in this example XXXX
plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:100,h_out$hi[order(h_out$hi[,1]),3],1:100,h_out$hi[order(h_out$hi[,1]),4])

## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
gc_out<-est_genocl(Gx=GlikHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="glik",ploidy="diploid",hier=TRUE,n_iters=4000)

## how variable is introgression among loci, lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but think this is mostly a good idea
sz_out<-sum2zero(hmc=gc_out$gencline_hmc,transform=TRUE,ci=0.90)

## plot genomic clines for the 51 loci, first without the sum-to-zero constraint
## then with it... these differ more for some data sets than others
gencline_plot(center=gc_out$center[,1],v=gc_out$gradient,pdf=FALSE)
gencline_plot(center=sz_out$center[,1],v=sz_out$gradient,pdf=FALSE)

## summarize loci with credible deviations from genome-averge gradients, here the focus is
## specifically on steep clines indicative of loci introgressing less than the average
which(sz_out$v[,2] > 1) ## index for loci with credibly steep clines
sum(sz_out$v[,2] > 1) ## number of loci with credibly steep clines

## last, lets look at interspecific ancestry for the same data set, this can
## be especially informative about the types of hybrids present
q_out<-est_Q(Gx=GlikHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="glik",ploidy="diploid")

## plot the results
tri_plot(hi=q_out$hi[,1],Q10=q_out$Q10[,1],pdf=FALSE,pch=19)
## note that some individuals appear to be likely backcrosses (close to the outer lines of the triangles)
## but the inidividals with intermediate hybrid indexes are clearly not F1s but rather late generation hybrids
```

# Citations

The general hierarchical Bayesian model used for Bayesidan genomic cline analysis was described here:

[Gompert Z, Buerkle CA (2011) Bayesian estimation of genomic clines. Molecular Ecology, 20:2111-2127.](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2011.05074.x)

The current set of models based on the log-logistic function and using HMC were first described here:

[Fierno TJ, Semenov G, Dopman EB, Taylor SA, Larson EL, Gompert Z (2023) Quantitative analyses of coupling in hybrid zones. Cold Spring Harb Perspect Biol, a041434.](https://cshperspectives.cshlp.org/content/early/2023/09/21/cshperspect.a041434)

The following paper (in the works) describes this specific R package:



