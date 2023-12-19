#' Function to estimate parental allele frequencies
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference parental allele frequencies from genetic data.
#' @param G0 genetic data for parental reference set 0 in the form of a matrix for known genotypes (rows = individuals, columns = loci) or a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype).
#' @param G1 genetic data for parental reference set 1 formatted as described for G0.
#' @param model for genetic data, either 'genotype' for known gentoypes or 'glik' for genotype likelihoods.
#' @param ploidy species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.
#' @param pldat list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid); matrix 2 in the list is for parental reference set 0 and matrix 3 is for parental reference set 1 (matriix 1 is reserved for the hybrids, which are not included in this analysis).
#' @param n_chains number of HMC chains for posterior inference.
#' @param n_iters A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.
#' @param p_warmup proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.
#' @param n_thin positive integer, save every n_thin HMC iterations for the posterior, default is 1.
#' @param n_cores number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if available).
#'
#' @details
#' Parental allele frequencies can be inferred based on known genotypes (model = 'genotype') or genotype likelihoods (model = 'glik'). Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, you can use 0 and 1 or 0 and 2. Genotype likelihoods should be on their natural scale (not phred scaled) and the values for each locus for an individual should sum to 1. The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, you must use the 0 and 2 matrixes to store the likelihoods of the two possible states.
#' @details
#' Ploidy data are only required for the mixed ploidy data. In this case, there should be a list of matrixes (the 1st matrix is for the hybrid so not used here), including one for each parent (2nd and 3rd matrixes, with parent 0 first). The matrixes indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1).
#'
#' @return A list of parameter estimates and full HMC results from stan. Parameter estimates are provided as a point estimate (median of the posterior) and 95% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. In this case, there are two matrixes, p0 and p1, for parent 0 and parent 1 allele frequencies. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default. 
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC 
#'
#' @export
est_p<-function(G0=NULL,G1=NULL,model="genotype",ploidy="diploid",pldat=NULL,
	hier=TRUE,SDc=NULL,SDv=NULL,n_chains=4,n_iters=2000,p_warmup=0.5,n_thin=1,n_cores=NULL){

        ## get or set number of cores for HMC
        if(is.null(n_cores)){
                options(mc.cores = parallel::detectCores())
        } else{
                mc.cores<-n_cores
        }

        ## determine number of warmup iterations
        n_warmup<-floor(p_warmup * n_iters)


	if(model=="genotype" & ploidy=="diploid"){
		dat<-list(L=dim(G0)[2],N0=dim(G0)[1],N1=dim(G1)[1],G0=G0,G1=G1)
		fit<-rstan::sampling(stanmodels$p,data=dat)
		p0<-t(apply(rstan::extract(fit,"P0")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		p1<-t(apply(rstan::extract(fit,"P1")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		pout<-list(p0=p0,p1=p1,p_hmc=fit)

		return(pout)
	} else if(model=="glik" & ploidy=="diploid"){
		dat<-list(L=dim(G0[[1]])[2],N0=dim(G0[[1]])[1],N1=dim(G1[[1]])[1],
		G00=G0[[1]],G10=G1[[1]],G01=G0[[2]],G11=G1[[2]],G02=G0[[3]],G12=G1[[3]])
		fit<-rstan::sampling(stanmodels$p_gl,data=dat)
		p0<-t(apply(rstan::extract(fit,"P0")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		p1<-t(apply(rstan::extract(fit,"P1")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		pout<-list(p0=p0,p1=p1,p_hmc=fit)
	
	} else if(model=="genotype" & ploidy=="mixed"){
		dat<-list(L=dim(G0)[2],N0=dim(G0)[1],N1=dim(G1)[1],G0=G0,G1=G1,
		ploidy0=pldat[[2]],ploidy1=pldat[[3]])
		fit<-rstan::sampling(stanmodels$p_mix,data=dat)
		p0<-t(apply(rstan::extract(fit,"P0")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		p1<-t(apply(rstan::extract(fit,"P1")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		pout<-list(p0=p0,p1=p1,p_hmc=fit)
	
	} else if(model=="genotype" & ploidy=="mixed"){
		dat<-list(L=dim(G0[[1]])[2],N0=dim(G0[[1]])[1],N1=dim(G1[[1]])[1],
		G00=G0[[1]],G10=G1[[1]],G01=G0[[2]],G11=G1[[2]],G02=G0[[3]],G12=G1[[3]],
		ploidy0=pldat[[2]],ploidy1=pldat[[3]])
		fit<-rstan::sampling(stanmodels$p_gl_mix,data=dat)
		p0<-t(apply(rstan::extract(fit,"P0")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		p1<-t(apply(rstan::extract(fit,"P1")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		pout<-list(p0=p0,p1=p1,p_hmc=fit)
	
	} else{
		stop("invalid model or ploidy specified")
	}
}
