# bgchm R package

This is an R package for Bayesian analyses of population genomic data from hybrid zones, including Bayesian genomic cline analysis, estimation of hybrid indexes and ancestry classes, some geographic cline analyses, and accessory plotting functions. This package using Hamiltonian Monte Carlo (HMC) for sampling posterior distributions, with HMC sampling implemented via [Stan](https://mc-stan.org).

# Installation

You can install this package within R directly from GitHub. Compiling the code requires a C++ compiler and associated components (see [Configuring C++ Toolchain](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#configuring-c-toolchain) from `rstan`).

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

## Fit genomic clines for an example data set with known genotypes

This data set comprises 51 ancestry-informative loci, 100 putative hybrids and 50 individuals represntative of each parental species. All loci are diploid. The data were simulated with [dfuse](https://cbuerkle.bitbucket.io/software/dfuse/) assuming 110 demes, m = 0.1 between neighboring demes, and 10 affecting hybrid fitness via underdominance (the underdominance model for `dfuse` is described in [Fierno et al. 2023](https://cshperspectives.cshlp.org/content/early/2023/09/21/cshperspect.a041434).

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
## sorted by hybrid index, just a nice way to visualize that in this example we have
## few hybrids with intermediate hybrid indexes

plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:100,h_out$hi[order(h_out$hi[,1]),3],1:100,h_out$hi[order(h_out$hi[,1]),4])

## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
gc_out<-est_genocl(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,n_iters=4000)

## how variable is introgression among loci? Lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but this is mostly a good idea
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

## Fit genomic clines for an example data set with genotype likelihoods

This data set comprises 51 ancestry-informative loci, 100 putative hybrids and 50 individuals represntative of each parental species. All loci are diploid. The data were simulated with `dfuse` assuming 110 demes, m = 0.1 between neighboring demes, and 100 loci weakly affecting hybrid fitness via underdominance.

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

## how variable is introgression among loci? Let's look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but this is mostly a good idea
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

## Fit genomic clines for an example data set with many loci

This example demonstrates how to fit clines for many loci in a scalable, parallelizable manner. The first step is to use a representative subset of loci (1000 or fewer) to estimate hybrid indexes and cline parameter standard deviations. Once this is done, we can analyze the remainder of the loci indpendently, which means this can be done in parallel. In this example, I will parallelize this part of the analysis with a bash script that fits sets of 10,000 SNPs in 20 jobs (200,000 SNPs total) where up to 10 jobs run at once. Each job uses 4 cores (for the four HMC chains) and thus my example assumes you have a computer with 10 x 4 = 40 cores. This can of course be adjusted to fit the details of your computational resources. Each set of 10,000 SNPs is saved in a unique rda object, and these are then loaded and combined in R in a final step.

### Step 1: Estimate the cline SDs and hybrid indexes
```R
## load the data set
data(manyloci)
## this includes three objects, G200kHybrids, G200kP0, and G200kP1

## select 1000 loci at random to estimate hybrid indexes and cline SDs
L<-dim(G200kHybrids)[2]
rset<-sample(1:L,1000,replace=FALSE)

## estimate hybrid indexes, uses default HMC settings
## and estimates allele frequencies on the fly
h_out<-est_hi(Gx=G200kHybrids[,rset],G0=G200kP0[,rset],G1=G200kP1[,rset],model="genotype",ploidy="diploid")
## let's write the hybrid index estimates to a text file so we can easily access them later
write.table(file="h_est.txt",h_out$hi,row.names=FALSE,quote=FALSE)

## fit a hierarchical genomic cline model for subset of loci using the estimated
## hybrid indexes, estimate parental allele frequencies on the fly
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
## our main goal here is to estimate the cline SDs
gc_out<-est_genocl(Gx=G200kHybrids[,rset],G0=G200kP0[,rset],G1=G200kP1[,rset],H=h_out$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,n_iters=4000)

## here are the cline SDs
## we will need these for the next step
gc_out$SDc
#      50%        5%       95%
#0.2027722 0.1853855 0.2202991
gc_out$SDv
#      50%        5%       95%
#0.2717172 0.2629907 0.2803867
```

### Step 2: Estimate clines for all of the loci in parallel

Here is the bash script, this could be run interactively or put in a SLURM script or equivalent to run on a cluster (that is what I normally do).
```bash
#!/bin/sh
## max = number of jobs to run at once
max=10
## total = total number of jobs
total=20
count=0

for ((j=1; j<=20; j++))
do
    Rscript --vanilla fitSnps.R "$j" &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
```
Here is the related R script, `fitSnps.R`, which is called by the bash script. I am assuming everything is in the same directory and thus not providing paths for any files.

```R
## load the data set
data(manyloci)
## this includes three objects, G200kHybrids, G200kP0, and G200kP1

## read in prior estimates of hybrid index
h<-read.table("h_est.txt",header=TRUE)

## enter point estimates of cline standard deviations inferred from the subset of loci
## these are the actual values I estimated, if you run this example yours will likely be
## similar but not identical
sdc<-0.20
sdv<-0.27

## selection the appropriate subset of SNPs
## using in this case 10,000 SNP batches and
## getting the batch number from the command line
## through the bash script above

## batch size
bsize<-10000

## batch number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

## bounds for this batch
## if the last batch is incomplete, you could set some maximum
## upper bound (ub) here
lb<-(k-1)*bsize + 1
ub<-lb+bsize-1

## number of loci and hybrids for this analysis
L<-length(lb:ub)
N<-length(h)

## subset genotype objects
sGhyb<-G200kHybrids[lb:ub,]
sGP0<-G200kP0[lb:ub,]
sGP1<-G200kP1[lb:ub,]

## in this example, I am just saving the center and gradient,
## point estimates and CIs, for each locus
## one snp at a time, L rows for loci, 3 columns for point est and CIs
onev<-matrix(NA,nrow=L,ncol=3)
onec<-matrix(NA,nrow=L,ncol=3)

## loop over loci, fitting one at a time
## using the pre-estimated values of sdc and sdv
## save cline parameter estimates for each locus
for(i in 1:L){
	out<-est_genocl(Gx=sGhyb[i,],G0=sGP0[i,],G1=sGP1[i,],H=h[,1],model="genotype",
			ploidy="diploid",hier=FALSE,SDc=sdc,SDv=sdv)
	onev[i,]<-out$gradient
	onec[i,]<-out$center
}

## save the R object with batch ID included
out<-paste("clinesOut",k,".rda",sep="")
save(list=ls(),file=out)
```
### Step 3: Combine estimates from each batch

In this last step, we will combine the results from each data subset (batch) in a single R object. Although not shown in detail here, we can then generate any plots, impose the sum-to-zero constraint, etc.

```R
## get a list of all of the rda files
cf<-list.files(pattern="clinesOut")

## load the first rda file and create our storage objects
load(cf[[1]])
Lt<-dim(G200kHybrids)[2] ## number of loci
## object for gradient = v and center = c
## one row per locus, 3 columns for point estimate and CIs
v_est<-matrix(NA,nrow=Lt,ncol=3)
c_est<-matrix(NA,nrow=Lt,ncol=3)

## using lb and ub stored in rda file
v_est[lb:ub,]<-onev
c_est[lb:ub,]<-onec

## loop over all of the rda files, adding the estimates
for(cfi in 2:length(cf)){
        load(cf[cfi])
        v_est[lb:ub,]<-onev
        c_est[lb:ub,]<-onec
}
## save object
save(list=ls(),file="combinedClines.rda")

## now impose s2z constraints, make plots, etc.
```

## Fit geographic clines for example data set of allele frequencies

This data set comprises allele frequencies for 51 diploid loci from 110 demes. The data were simulated with m = 0.1 between neighboring demes and 10 loci affecting hybrid fitness via underdominance.
```R
## load the data set
data(pfreqs)
## this includes one object, a matrix P with allele frequencies
## 110 rows = demes, 51 columns = loci

## use standardized deme numbers as geographic coordinates
x<-1:110
geo<-(x-mean(x))/sd(x)

## fit the geographic cline model
o<-est_geocl(P=P,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE,n_iters=5000)

## plot clines on logit scale, which should be linear
plot(geo,o$cent[1,1] + o$slope[1,1] * geo,type='l',ylim=c(-15,15),ylab="Logit allele frequency",xlab="Deme number",
	axes=FALSE)
axis(1,at=geo[seq(5,110,5)],x[seq(5,110,5)])
axis(2)
box()
for(i in 2:51){
	lines(geo,o$cent[i,1] + o$slope[i,1] * geo)
}

```

# Citations

The general hierarchical Bayesian model used for Bayesidan genomic cline analysis was described here:

[Gompert Z, Buerkle CA (2011) Bayesian estimation of genomic clines. Molecular Ecology, 20:2111-2127.](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2011.05074.x)

The current set of models based on the log-logistic function and using HMC were first described here:

[Fierno TJ, Semenov G, Dopman EB, Taylor SA, Larson EL, Gompert Z (2023) Quantitative analyses of coupling in hybrid zones. Cold Spring Harb Perspect Biol, a041434.](https://cshperspectives.cshlp.org/content/early/2023/09/21/cshperspect.a041434)

The following paper (in the works) describes this specific R package:
