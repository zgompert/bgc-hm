% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_p.R
\name{est_p}
\alias{est_p}
\title{Function to estimate parental allele frequencies}
\usage{
est_p(
  G0 = NULL,
  G1 = NULL,
  model = "genotype",
  ploidy = "diploid",
  pldat = NULL,
  n_chains = 4,
  n_iters = 2000,
  p_warmup = 0.5,
  n_thin = 1,
  n_cores = NULL,
  HMC = FALSE,
  full = TRUE
)
}
\arguments{
\item{G0}{genetic data for parental reference set 0 in the form of a matrix for known genotypes (rows = individuals, columns = loci) or a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype).}

\item{G1}{genetic data for parental reference set 1 formatted as described for G0.}

\item{model}{for genetic data, either 'genotype' for known gentoypes or 'glik' for genotype likelihoods.}

\item{ploidy}{species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.}

\item{pldat}{list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid); matrix 2 in the list is for parental reference set 0 and matrix 3 is for parental reference set 1 (matrix 1 is reserved for the hybrids, which are not included in this analysis).}

\item{n_chains}{number of HMC chains for posterior inference.}

\item{n_iters}{A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.}

\item{p_warmup}{proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.}

\item{n_thin}{positive integer, save every n_thin HMC iterations for the posterior, default is 1.}

\item{n_cores}{number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if available).}

\item{HMC}{boolean denotes whether to use HMC (TRUE) or an analytical solution (FALSE) for the posterior for the genotype likelihood model (the analytical solution is always used for known genotypes), default is TRUE.}

\item{full}{boolean denoting whether (TRUE, the default) or not (FALSE) to return the HMC Stan object with the full set of samples from the posterior (only relevant if HMC = TRUE).}
}
\value{
A list of parameter estimates and full HMC results from stan (only when HMC is used). Parameter estimates are provided as a point estimate (median of the posterior), 90\% equal-tail probability intervals (5th and 95th quantiles of the posterior distribution), and 95\% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. In this case, there are two matrixes, p0 and p1, for parent 0 and parent 1 allele frequencies. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default.
}
\description{
Bayesian inference parental allele frequencies from genetic data.
}
\details{
Parental allele frequencies can be inferred based on known genotypes (model = 'genotype') or genotype likelihoods (model = 'glik'). With known genotypes, the posterior is computed exactly (it is a beta distribution given the binomial likelihood and a beta prior). With genotype likelihoods, the user can decide to compute the exact posterior based on a point estimate of the allele counts from the genotype likelihoods (HMC=FALSE), or to conduct the full Hamiltonian Monte Carlo analyses (HMC=TRUE). The latter provides a better characterization of uncertainty, but will take much more time. Thus, the analytical approach should be chosen for large (more than a few thousand loci) data sets.

Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, use 0 and 1. Genotype likelihoods should be on their natural scale (not phred scaled) and the values for each locus for an individual should sum to 1. The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, use the 0 and 1 matrixes to store the likelihoods of the two possible states. For all models, missing data can be encoded by setting the ploidy for an individual/locus to 0 (this indicates no information, whereas genotype likelihoods encode uncertainty in genotypes) and the genotype to NA.

Ploidy data are only required for the mixed ploidy data. In this case, there should be a list of matrixes (the 1st matrix is for the hybrid so not used here), including one for each parent (2nd and 3rd matrixes, with parent 0 first). The matrixes indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1) (use 0 for missing data).
}
\examples{
\dontrun{
## load the data set
data(genotypes)
## this includes three objects, GenHybrids, GenP0, and GenP1

## estimate parental allele frequencies, uses the analytical solution to the posterior
p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")
}
}
\seealso{
'rstan::stan' for details on HMC with stan and the rstan HMC
}
