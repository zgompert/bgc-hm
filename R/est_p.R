#' Function to estimate parental allele frequencies
#'
#' Bayesian inference parental allele frequencies from genetic data.
#' @param G0 genetic data for parental reference set 0 in the form of a matrix for known genotypes (rows = individuals, columns = loci) or a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype).
#' @param G1 genetic data for parental reference set 1 formatted as described for G0.
#' @param model for genetic data, either 'genotype' for known gentoypes or 'glik' for genotype likelihoods.
#' @param ploidy species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.
#' @param pldat list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid); matrix 2 in the list is for parental reference set 0 and matrix 3 is for parental reference set 1 (matrix 1 is reserved for the hybrids, which are not included in this analysis).
#' @param n_chains number of HMC chains for posterior inference.
#' @param n_iters A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.
#' @param p_warmup proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.
#' @param n_thin positive integer, save every n_thin HMC iterations for the posterior, default is 1.
#' @param n_cores number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if available).
#' @param HMC boolean denotes whether to use HMC (TRUE) or an analytical solution (FALSE) for the posterior for the genotype likelihood model (the analytical solution is always used for known genotypes), default is TRUE.
#' @param full boolean denoting whether (TRUE, the default) or not (FALSE) to return the HMC Stan object with the full set of samples from the posterior (only relevant if HMC = TRUE).
#'
#' @details
#' Parental allele frequencies can be inferred based on known genotypes (model = 'genotype') or genotype likelihoods (model = 'glik'). With known genotypes, the posterior is computed exactly (it is a beta distribution given the binomial likelihood and a beta prior). With genotype likelihoods, the user can decide to compute the exact posterior based on a point estimate of the allele counts from the genotype likelihoods (HMC=FALSE), or to conduct the full Hamiltonian Monte Carlo analyses (HMC=TRUE). The latter provides a better characterization of uncertainty, but will take much more time. Thus, the analytical approach should be chosen for large (more than a few thousand loci) data sets. 
#' @details
#' Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, use 0 and 1. Genotype likelihoods should be on their natural scale (not phred scaled) and the values for each locus for an individual should sum to 1. The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, use the 0 and 1 matrixes to store the likelihoods of the two possible states. For all models, missing data can be encoded by setting the ploidy for an individual/locus to 0 (this indicates no information, whereas genotype likelihoods encode uncertainty in genotypes).
#' @details
#' Ploidy data are only required for the mixed ploidy data. In this case, there should be a list of matrixes (the 1st matrix is for the hybrid so not used here), including one for each parent (2nd and 3rd matrixes, with parent 0 first). The matrixes indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1) (use 0 for missing data).
#'
#' @return A list of parameter estimates and full HMC results from stan (only when HMC is used). Parameter estimates are provided as a point estimate (median of the posterior), 90% equal-tail probability intervals (5th and 95th quantiles of the posterior distribution), and 95% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. In this case, there are two matrixes, p0 and p1, for parent 0 and parent 1 allele frequencies. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default. 
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC 
#'
#' @export
#' @examples
#'\dontrun{
#' ## load the data set
#' data(genotypes)
#' ## this includes three objects, GenHybrids, GenP0, and GenP1
#'
#' ## estimate parental allele frequencies, uses the analytical solution to the posterior
#' p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")
#' }
est_p<-function(G0=NULL,G1=NULL,model="genotype",ploidy="diploid",pldat=NULL,
	n_chains=4,n_iters=2000,p_warmup=0.5,n_thin=1,n_cores=NULL,HMC=FALSE,
	full=TRUE){

        ## get or set number of cores for HMC
        if(is.null(n_cores)){
                options(mc.cores = parallel::detectCores())
        } else{
                mc.cores<-n_cores
        }

        ## determine number of warmup iterations
        n_warmup<-floor(p_warmup * n_iters)
	
	## jeffery prior for allele frequencies (alpha0 and beta0 for beta prior)
	a0<-.5
	b0<-.5

	mqbeta<-function(p=NA,shape1=NA,shape2=NA){
		o<-cbind(qbeta(p[1],shape1=a,shape2=b),
		      qbeta(p[2],shape1=a,shape2=b),
		      qbeta(p[3],shape1=a,shape2=b),
		      qbeta(p[4],shape1=a,shape2=b),
		      qbeta(p[5],shape1=a,shape2=b))
		colnames(o)<- paste(p*100,"%",sep="")
		return(o)
	}

	if(model=="genotype" & ploidy=="diploid"){
		G0<-as.matrix(G0)
		G1<-as.matrix(G1)
		## solve for posterior
		y<-apply(G0,2,sum)
		n<-rep(dim(G0)[1] *2, dim(G0)[2])
		a<-y+a0
		b<-n-y+b0
		p0<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
		y<-apply(G1,2,sum)
		n<-rep(dim(G1)[1] *2, dim(G1)[2])
		a<-y+a0
		b<-n-y+b0
		p1<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
		pout<-list(p0=p0,p1=p1)

	} else if(model=="glik" & ploidy=="diploid"){
		for(k in 1:3){
			G0[[k]]<-as.matrix(G0[[k]])
			G1[[k]]<-as.matrix(G1[[k]])
		}
		if(HMC==TRUE){ ## use HMC
			dat<-list(L=dim(G0[[1]])[2],N0=dim(G0[[1]])[1],N1=dim(G1[[1]])[1],
			GL00=G0[[1]],GL10=G1[[1]],GL01=G0[[2]],GL11=G1[[2]],GL02=G0[[3]],GL12=G1[[3]])
			fit<-rstan::sampling(stanmodels$p_gl,data=dat,
				iter=n_iters,warmup=n_warmup,thin=n_thin)
			p0<-t(apply(rstan::extract(fit,"P0")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
			p1<-t(apply(rstan::extract(fit,"P1")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
			## create a list with parameter estimates plus full hmc object
			if(full==TRUE){
				pout<-list(p0=p0,p1=p1,p_hmc=fit)
			} else{
				pout<-list(p0=p0,p1=p1)
			}
		} else{ ## use analytical
			## solve for posterior
			y<-apply(G0[[2]]+G0[[3]]*2,2,sum)
			n<-rep(dim(G0[[1]])[1] *2, dim(G0[[1]])[2])
			a<-y+a0
			b<-n-y+b0
			p0<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
			y<-apply(G1[[2]]+G1[[3]]*2,2,sum)
			n<-rep(dim(G1[[1]])[1] *2, dim(G1[[1]])[2])
			a<-y+a0
			b<-n-y+b0
			p1<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
			pout<-list(p0=p0,p1=p1)
		}
	} else if(model=="genotype" & ploidy=="mixed"){
		G0<-as.matrix(G0)
		G1<-as.matrix(G1)
		pldat[[2]]<-as.matrix(pldat[[2]])
		pldat[[3]]<-as.matrix(pldat[[3]])
		## solve for posterior
		y<-apply(G0,2,sum,na.rm=TRUE)
		n<-apply(pldat[[2]],2,sum)
		a<-y+a0
		b<-n-y+b0
		p0<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
		y<-apply(G1,2,sum,na.rm=TRUE)
		n<-apply(pldat[[3]],2,sum)
		a<-y+a0
		b<-n-y+b0
		p1<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
		pout<-list(p0=p0,p1=p1)
	
	} else if(model=="glik" & ploidy=="mixed"){
		for(k in 1:3){
			G0[[k]]<-as.matrix(G0[[k]])
			G1[[k]]<-as.matrix(G1[[k]])
		}
		pldat[[2]]<-as.matrix(pldat[[2]])
		pldat[[3]]<-as.matrix(pldat[[3]])
		
		if(HMC==TRUE){ ## use HMC	
			dat<-list(L=dim(G0[[1]])[2],N0=dim(G0[[1]])[1],N1=dim(G1[[1]])[1],
			GL00=G0[[1]],GL10=G1[[1]],GL01=G0[[2]],GL11=G1[[2]],GL02=G0[[3]],GL12=G1[[3]],
			ploidy0=pldat[[2]],ploidy1=pldat[[3]])
			fit<-rstan::sampling(stanmodels$p_gl_mix,data=dat,
				iter=n_iters,warmup=n_warmup,thin=n_thin)
			p0<-t(apply(rstan::extract(fit,"P0")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
			p1<-t(apply(rstan::extract(fit,"P1")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
			## create a list with parameter estimates plus full hmc object
			if(full==TRUE){
				pout<-list(p0=p0,p1=p1,p_hmc=fit)
			} else{
				pout<-list(p0=p0,p1=p1)
			}
		} else { ## use analytical
			## solve for posterior
			y<-apply(G0[[2]]+G0[[3]]*2,2,sum,na.rm=TRUE)
			n<-apply(pldat[[2]],2,sum)
			a<-y+a0
			b<-n-y+b0
			p0<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
			y<-apply(G1[[2]]+G1[[3]]*2,2,sum,na.rm=TRUE)
			n<-n<-apply(pldat[[3]],2,sum)
			a<-y+a0
			b<-n-y+b0
			p1<-mqbeta(c(.5,.025,.05,.95,.975),shape1=a,shape2=b)
			pout<-list(p0=p0,p1=p1)
		}
	
	} else{
		stop("invalid model or ploidy specified")
	}
	return(pout)
}
