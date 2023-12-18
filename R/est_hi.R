#' Function to estimate hybrid index
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of hybrid indexes from genetic data.
#' @param Gx genetic data for putative hybrids in the form of a matrix for known genotypes (rows = individuals, columns = loci), a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype), or  matrix of ancestry for the ancestry model (rows = individuals, columns = loci).
#' @param G0 genetic data for parental reference set 0 formatted as described for Gx.
#' @param G1 genetic data for parental reference set 1 formatted as described for Gx.
#' @param p0 vector of allele frequencies for parental reference set 0 (one entry per locus).
#' @param p1 vector allele frequencies for parental reference set 1 (one entry per locus).
#' @param H vector of hybrid indexes for the putative hybrids (one entry per locus).
#' @param model for genetic data, either 'genotype' for known gentoypes, 'glik' for genotype likelihoods, or 'ancestry' for known ancestry.
#' @param ploidy species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.
#' @param pldat matrix or list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid).
#' @param n_chains number of HMC chains for posterior inference.
#' @param n_iters A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.
#' @param p_warmup proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.
#' @param n_thin positive integer, save every n_thin HMC iterations for the posterior, default is 1.
#' @param n_cores number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if available).
#'
#' @return A list of parameter estimates (hi = hybrid indexes) and full HMC results from stan. Parameter estimates are provided as a point estimate (median of the posterior) and 95% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default.
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC output object.
#'
#' @references
#' Gompert Z, et al. 2024. Bayesian hybrid zone analyses with Hamiltonian Monte Carlo in R. Manuscript in preparation
#' @export
est_hi<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,model="genotype",ploidy="diploid",pldat=NULL,
		 n_chains=4,n_iters=2000,p_warmup=0.5,n_thin=1,n_cores=NULL){

        ## get or set number of cores for HMC
        if(is.null(n_cores)){
                options(mc.cores = parallel::detectCores())
        } else{
                mc.cores<-n_cores
        }

        ## determine number of warmup iterations
        n_warmup<-floor(p_warmup * n_iters)

	if(is.list(Gx)==TRUE & model!="glik"){
                stop("List input is only valide for the glik model")
        } else if(is.list(Gx)==FALSE & model=="glik"){
                stop("List input required for the genotype likelihood model")
        } else if(! model %in% c("genotype","glik","ancestry")){
            stop("Unknown model, must be genotype, glik, or ancestry")
        } else if(! ploidy %in% c("diploid","mixed")){
            stop("Only diploid and mixed are accepted for the ploidy argument")
        } else if(ploidy=="mixed" & is.null(pldat)){
            stop("pldat cannot be NULL for mixed ploitdy")
        }

        ## estimate parental allele frequencies if not provided
        if(is.null(p0) | is.null(p1)){
                po<-est_p(G0=G0,G1=G1,model=model,ploidy=ploidy,pldat=pldat)
                p0<-po$p0[,1]
                p1<-po$p1[,1]
        }

        ## now have population allele frequencies, just need hybrid ploidy data if it
        ## was provided as a list
        if(is.list(pldat)==TRUE){
                pldat<-pldat[[1]]
        }

	## fit the appropriate hybrid index model based on options
	if(model=="genotype" & ploidy=="diploid"){
		## diploid with known genotypes
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1)
		fit<-rstan::sampling(stanmodels$hi,data=dat,
                   iter=n_iters,warmup=p_warmup,thin=n_thin)
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		hiout<-list(hi=hi,hi_hmc=fit)
		return(hiout)
	} else if(model=="glik" & ploidy=="diploid"){
		## diploid with genotype likelihoods
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]]
			  ,P0=p0,P1=p1)
		fit<-rstan::sampling(stanmodels$hi_gl,data=dat,
                   iter=n_iters,warmup=p_warmup,thin=n_thin)
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		hiout<-list(hi=hi,hi_hmc=fit)
		return(hiout)
	} else if(model=="ancestry" & ploidy=="diploid"){
		## diploid with known ancestry
		return(hiout)
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx)
		fit<-rstan::sampling(stanmodels$hi_z,data=dat,
                   iter=n_iters,warmup=p_warmup,thin=n_thin)
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		hiout<-list(hi=hi,hi_hmc=fit)
	} else if(model=="genotype" & ploidy=="mixed"){
		## diploid with known genotypes
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1,ploidy=pldat)
		fit<-rstan::sampling(stanmodels$hi_mix,data=dat,
                   iter=n_iters,warmup=p_warmup,thin=n_thin)
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		## create a list with parameter estimates plus full hmc object
		hiout<-list(hi=hi,hi_hmc=fit)
		return(hiout)
	} ## need to add mixed glik and mixed ancestry still
	
}
