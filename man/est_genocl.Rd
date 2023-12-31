% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_gencline.R
\name{est_genocl}
\alias{est_genocl}
\title{Estimate genomic clines using Bayesian HMC}
\usage{
est_genocl(
  Gx = NULL,
  G0 = NULL,
  G1 = NULL,
  p0 = NULL,
  p1 = NULL,
  H = NULL,
  model = "genotype",
  ploidy = "diploid",
  pldat = NULL,
  hier = TRUE,
  SDc = NULL,
  SDv = NULL,
  n_chains = 4,
  n_iters = 2000,
  p_warmup = 0.5,
  n_thin = 1,
  n_cores = NULL
)
}
\arguments{
\item{Gx}{genetic data for putative hybrids in the form of a matrix for known genotypes (rows = individuals, columns = loci), a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype), or  matrix of ancestry for the ancestry model (rows = individuals, columns = loci).}

\item{G0}{genetic data for parental reference set 0 formatted as described for Gx.}

\item{G1}{genetic data for parental reference set 1 formatted as described for Gx.}

\item{p0}{vector of allele frequencies for parental reference set 0 (one entry per locus).}

\item{p1}{vector allele frequencies for parental reference set 1 (one entry per locus).}

\item{H}{vector of hybrid indexes for the putative hybrids (one entry per locus).}

\item{model}{for genetic data, either 'genotype' for known gentoypes, 'glik' for genotype likelihoods, or 'ancestry' for known ancestry.}

\item{ploidy}{species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.}

\item{pldat}{matrix or list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid).}

\item{hier}{Boolean, fit hierarchical model (TRUE) that estimates cline SDs or non-hierarchical model that assumes cline SDs are known (FALSE).}

\item{SDc}{known cline center SD on logit scale.}

\item{SDv}{known cline gradient SD on log10 scale.}

\item{n_chains}{number of HMC chains for posterior inference.}

\item{n_iters}{A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.}

\item{p_warmup}{proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.}

\item{n_thin}{positive integer, save every n_thin HMC iterations for the posterior, default is 1.}

\item{n_cores}{number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if available).}
}
\value{
A list of parameter estimates and full HMC results from stan, this includes cline parameters (center and gradient), and, for hierarchical models, standard deviations describing variability in clines across loci (SDc and SDv). Parameter estimates are provided as a point estimate (median of the posterior) and 95\% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default.
}
\description{
Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of genomic clines from genetic data. This fits a hierarchical log-logistic model for genomic clines with two key parameters, cline center and cline gradient (i.e., slope, inversely proportional to cline width). The model also estimates the cline standard deviations (SDs) across loci, SDc = variation in logit centers and SDv = variation in log10 gradients.
}
\details{
Clines can be inferred based on known genotypes (model = 'genotype'), genotype likelihoods, (model = 'glik') or known (estimated) ancestry (model = 'ancestry'). Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, you can use 0 and 1 or 0 and 2. Genotype likelihoods should be on their natural scale (not phred scaled) and the values for each locus for an individual should sum to 1. The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, you must use the 0 and 2 matrixes to store the likelihoods of the two possible states. For the ancestry model, clines are inferred directly from known local (locus-specific) ancestry rather than from genotype data. Users are free to use whatever software they prefer for local ancestry inference (many exist). In this case, each entry in the individual (rows) by locus (columns) matrix should denote the number of gene copies inherited from parental population 1 (where pure parent 1 corresponds with a hybrid index of 1 and pure parent 0 corresponds with a hybrid index of 0). Haploids can be encoded using 0 and 1 or 0 and 2 (this is treated equivalently).

Hybrid indexes must be provided. Hybrid genetic (or ancestry) data are also always required. For genotype or genotype likelihood models, users must either provide pre-estimated parental allele frequencies or parent genetic (genotypes or genotype likelihoods) that can be used to infer allele frequencies. Parental data are not required for the ancestry model.

Ploidy data are only required for the mixed ploidy data. In this case, there should be one matrix for the hybrids or a list of matrixes for the hybrids (1st matrix) and each parent (2nd and 3rd matrixes, with parent 0 first). The latter is required for the genotype or genotype likelihood models if parental allele frequencies are not provided. The matrixes indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1).
}
\references{
Gompert Z, et al. 2024. Bayesian hybrid zone analyses with Hamiltonian Monte Carlo in R. Manuscript in preparation
}
\seealso{
'rstan::stan' for details on HMC with stan and the rstan HMC output object.
}
