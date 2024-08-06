#' Estimate genomic clines using Bayesian HMC
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of genomic clines from genetic data. This fits a hierarchical logit-logistic model for genomic clines with two key parameters, cline center and cline gradient (i.e., slope, inversely proportional to cline width). The model also estimates the cline standard deviations (SDs) across loci, SDc = variation in logit centers and SDv = variation in log10 gradients. 
#' @param Gx genetic data for putative hybrids in the form of a matrix for known genotypes (rows = individuals, columns = loci), a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype), or  matrix of ancestry for the ancestry model (rows = individuals, columns = loci).
#' @param G0 genetic data for parental reference set 0 formatted as described for Gx.
#' @param G1 genetic data for parental reference set 1 formatted as described for Gx.
#' @param p0 vector of allele frequencies for parental reference set 0 (one entry per locus).
#' @param p1 vector allele frequencies for parental reference set 1 (one entry per locus).
#' @param H vector of hybrid indexes for the putative hybrids (one entry per locus).
#' @param model for genetic data, either 'genotype' for known gentoypes, 'glik' for genotype likelihoods, or 'ancestry' for known ancestry.
#' @param ploidy species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.
#' @param pldat matrix or list of matrixes with ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid).
#' @param hier Boolean, fit hierarchical model (TRUE) that estimates cline SDs or non-hierarchical model that assumes cline SDs are known (FALSE).
#' @param SDc known cline center SD on logit scale.
#' @param SDv known cline gradient SD on log10 scale.
#' @param estMu boolean, estimate genomic cline means (TRUE) or assume a prior mean of 0 (FALSE, the default). 
#' @param sd0 standard deviation for the normal prior on the cline standard deviations, default is 1.
#' @param mu0 standard deviation for the normal prior on the cline mean, default is 1.
#' @param n_chains number of HMC chains for posterior inference.
#' @param n_iters A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.
#' @param p_warmup proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.
#' @param n_thin positive integer, save every n_thin HMC iterations for the posterior, default is 1.
#' @param n_cores number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if more cores are available).
#'
#' @details
#' Clines can be inferred based on known genotypes (model = 'genotype'), genotype likelihoods, (model = 'glik') or known (estimated) ancestry (model = 'ancestry'). Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, use 0 and 1. Genotype likelihoods should be on their natural scale (not on the phred scaled) and the values for each locus for an individual should sum to 1 (i.e., the likelihoods are scaled to be probabilities). The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, you must use the 0 and 1 matrixes to store the likelihoods of the two possible states. For the ancestry model, clines are inferred directly from known local (locus-specific) ancestry rather than from genotype data. Users are free to use whatever software they prefer for local ancestry inference (many exist). In this case, each entry in the individual (rows) by locus (columns) matrix should denote the number of gene copies inherited from parental population 1 (where pure parent 1 corresponds with a hybrid index of 1 and pure parent 0 corresponds with a hybrid index of 0). Haploids can be encoded using 0 and 1. For all models, missing data can be encoded by setting the ploidy for an individual/locus to 0 (this indicates no information, whereas genotype likelihoods encode uncertainty in genotypes). 
#' @details
#' Hybrid indexes must be provided. Genetic (or ancestry) data for putative hybrids are also always required. For genotype or genotype likelihood models, users must either provide pre-estimated parental allele frequencies or parent genetic (genotypes or genotype likelihoods) that can be used to infer allele frequencies with est_p. Parental data are not required for the ancestry model.
#' @details
#' Ploidy data are only required for the mixed ploidy data. In this case, there should be one matrix for the hybrids or a list of matrixes for the hybrids (1st matrix) and each parent (2nd and 3rd matrixes, with parent 0 first). The latter is required for the genotype or genotype likelihood models if parental allele frequencies are not provided. The matrixes indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1) (use 0 for missing data).
#'
#' @return A list of parameter estimates and full HMC results from stan, this includes cline parameters (center and gradient), and, for hierarchical models, standard deviations describing variability in clines across loci (SDc and SDv). Mean values for cline center (Muc) and slope (Muv) are also provided if estMu is set to TRUE. Parameter estimates are provided as a point estimate (median of the posterior) and 90% equal-tail probability intervals (5th and 95th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default. 
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC output object.
#'
#' @references
#' Gompert Z, DeRaad D, Buerkle CA. A next generation of hierarchical Bayesian analyses of hybrid zones enables direct quantification of variation in introgression in R. bioRxiv 2024.03.29.587395.
#' @export
#' @examples
#'\dontrun{
#' ## load the data set
#' data(genotypes)
#' ## this includes three objects, GenHybrids, GenP0, and GenP1
#'
#' ## estimate parental allele frequencies
#' p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")
#'
#' ## estimate hybrid indexes, uses default HMC settings
#' ## and uses point estimates (posterior medians) of allele frequencies
#' h_out<-est_hi(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")
#'
#' ## fit a hierarchical genomic cline model for all 51 loci using the estimated
#' ## hybrid indexes and parental allele frequencies (point estimates)
#' ## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
#' gc_out<-est_genocl(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,n_iters=4000)
#' }
est_genocl<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,H=NULL,model="genotype",ploidy="diploid",pldat=NULL,
	hier=TRUE,SDc=NULL,SDv=NULL,estMu=FALSE,sd0=1,mu0=1,n_chains=4,n_iters=2000,p_warmup=0.5,n_thin=1,n_cores=NULL){
	
	## get or set number of cores for HMC
	if(is.null(n_cores)){
		options(mc.cores = parallel::detectCores())
	} else{
		mc.cores<-n_cores
	}
	
	## determine number of warmup iterations
	n_warmup<-floor(p_warmup * n_iters)
	
	## sanity tests to make sure options match up
	if(hier==FALSE & (is.null(SDc) | is.null(SDv))){
		stop("SDc and SDV required for hier==FALSE")
	} else if(is.list(Gx)==TRUE & model!="glik"){
		stop("List input is only valide for the glik model")
	} else if(is.list(Gx)==FALSE & model=="glik"){
		stop("List input required for the genotype likelihood model")
	} else if(! model %in% c("genotype","glik","ancestry")){
	    stop("Unknown model, must be genotype, glik, or ancestry") 
	} else if(! ploidy %in% c("diploid","mixed")){
	    stop("Only diploid and mixed are accepted for the ploidy argument") 
	} else if(ploidy=="mixed" & is.null(pldat)){
	    stop("pldat cannot be NULL for mixed ploitdy")
	} else if(estMu==TRUE & hier==FALSE){
		stop("genomic cline means can only be estimated using the hierarchical model")
	}
	
	## estimate parental allele frequencies if not provided
	if(is.null(p0) | is.null(p1)){
		po<-est_p(G0=G0,G1=G1,model=model,ploidy=ploidy,pldat=pldat,HMC=FALSE)
		p0<-po$p0[,1]
		p1<-po$p1[,1]
	}
	
	if(sum(p0==0&p1==0)+sum(p0==1&p1==1) > 0){
		message("one or more loci are fixed for the same allele in both parental source populations, this can cause problems\n")
	}
	
	## now have population allele frequencies, just need hybrid ploidy data if it
	## was provided as a list
	if(is.list(pldat)==TRUE){
		pldat<-pldat[[1]]
	}
	
	## fit appropriate genomic cline model based on options
	if(hier==TRUE & model=="genotype" & ploidy=="diploid" & estMu==FALSE){
	## hierachical model, known genotypes, diploids
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,sd0=sd0,
		init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline,data=dat,
		iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="genotype" & ploidy=="diploid" & estMu==TRUE){
	## hierachical model, known genotypes, diploids, estimate means
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,sd0=sd0,mu0=mu0,
		init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_mu,data=dat,
		iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		muc<-quantile(rstan::extract(fit,"muc")[[1]],probs=c(.5,.05,.95))
		muv<-quantile(rstan::extract(fit,"muv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,Muc=muc,Muv=muv,gencline_hmc=fit)
	} else if(hier==FALSE & model=="genotype" & ploidy=="diploid"){
	## known SD model, known genotypes, diploids
		if(is.matrix(Gx)){
			## generate initial values of cline parameters
			initf<-function(L=dim(Gx)[2],chain_id=1){
        			list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
			dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_sdk,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
			cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))

		} else{ ## one snp
			## generate initial values of cline parameters
			initf<-function(chain_id=1){
        			list(center=runif(1,.3,.7),v=runif(1,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))		
			dat<-list(L=1,N=length(Gx),G=Gx,H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_one,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)

	} else if(hier==TRUE & model=="glik" & ploidy=="diploid" & estMu==FALSE){
	## hierachical model, genotype likelihoods, diploids
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx[[1]])[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
		H=H,P0=p0,P1=p1,sd0=sd0,init=init_ll)

		if(dim(Gx[[1]])[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_gl,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="glik" & ploidy=="diploid" & estMu==TRUE){
	## hierachical model, genotype likelihoods, diploids, estimate means
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx[[1]])[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
		H=H,P0=p0,P1=p1,sd0=sd0,mu0=mu0,init=init_ll)

		if(dim(Gx[[1]])[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_gl_mu,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		mudc<-quantile(rstan::extract(fit,"muc")[[1]],probs=c(.5,.05,.95))
		mudv<-quantile(rstan::extract(fit,"muv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,Muc=muc,Muv=muv,gencline_hmc=fit)
	} else if(hier==FALSE & model=="glik" & ploidy=="diploid"){
	## known SD model, genotype likelihoods, diploids
		if(is.matrix(Gx[[1]])){
			## generate initial values of cline parameters
			initf<-function(L=dim(Gx[[1]])[2],chain_id=1){
        			list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
			dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],
			GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_sdk_gl,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
			cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))

		} else{ ## one snp
			## generate initial values of cline parameters
			initf<-function(chain_id=1){
        			list(center=runif(1,.3,.7),v=runif(1,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))		
			dat<-list(L=1,N=length(Gx[[1]]),GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
			H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_one_gl,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="ancestry" & ploidy=="diploid" & estMu==FALSE){
	## hierachical model, ancestry, diploids
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,sd0=sd0,init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_z,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="ancestry" & ploidy=="diploid" & estMu==TRUE){
	## hierachical model, ancestry, diploids, estimate means
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,sd0=sd0,mu0=mu0,init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_z_mu,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		muc<-quantile(rstan::extract(fit,"muc")[[1]],probs=c(.5,.05,.95))
		muv<-quantile(rstan::extract(fit,"muv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,Muc=muc,Muv=muv,gencline_hmc=fit)
	} else if(hier==FALSE & model=="ancestry" & ploidy=="diploid"){
	## known SD model, ancestry, diploids
		if(is.matrix(Gx)){
			## generate initial values of cline parameters
			initf<-function(L=dim(Gx)[2],chain_id=1){
        			list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
			dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,sc=SDc,sv=SDv,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_sdk_z,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
			cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))

		} else{ ## one snp
			## generate initial values of cline parameters
			initf<-function(chain_id=1){
        			list(center=runif(1,.3,.7),v=runif(1,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))		
			dat<-list(L=1,N=length(Gx),Z=Gx,H=H,sc=SDc,sv=SDv,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_one_z,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)

	} else if(hier==TRUE & model=="genotype" & ploidy=="mixed" & estMu==FALSE){
	## hierachical model, known genotypes, mixed ploidy
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,ploidy=pldat,sd0=sd0,
		init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="genotype" & ploidy=="mixed" & estMu==TRUE){
	## hierachical model, known genotypes, mixed ploidy, estimate means
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,ploidy=pldat,sd0=sd0,mu0=mu0,
		init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_mix_mu,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		muc<-quantile(rstan::extract(fit,"muc")[[1]],probs=c(.5,.05,.95))
		muv<-quantile(rstan::extract(fit,"muv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,Muc=muc,Muv=muv,gencline_hmc=fit)
	} else if(hier==FALSE & model=="genotype" & ploidy=="mixed"){
	## known SD model, known genotypes, mixed ploidy
		if(is.matrix(Gx)){
			## generate initial values of cline parameters
			initf<-function(L=dim(Gx)[2],chain_id=1){
        			list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
			dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,
			ploidy=pldat,init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_sdk_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
			cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))

		} else{ ## one snp
			## generate initial values of cline parameters
			initf<-function(chain_id=1){
        			list(center=runif(1,.3,.7),v=runif(1,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))		
			dat<-list(L=1,N=length(Gx),G=Gx,H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,ploidy=pldat,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_one_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="glik" & ploidy=="mixed" & estMu==FALSE){
	## hierachical model, genotype likelihoods, mixed ploidy
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
		H=H,P0=p0,P1=p1,ploidy=pldat,sd0=sd0,init=init_ll)

		if(dim(Gx[[1]])[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_gl_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="glik" & ploidy=="mixed" & estMu==TRUE){
	## hierachical model, genotype likelihoods, mixed ploidy, estimate means
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
		H=H,P0=p0,P1=p1,ploidy=pldat,sd0=sd0,mu0=mu0,init=init_ll)

		if(dim(Gx[[1]])[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_gl_mix_mu,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		muc<-quantile(rstan::extract(fit,"muc")[[1]],probs=c(.5,.05,.95))
		muv<-quantile(rstan::extract(fit,"muv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,Muc=muc,Muv=muv,gencline_hmc=fit)
	} else if(hier==FALSE & model=="glik" & ploidy=="mixed"){
	## known SD model, genotype likelihoods, mixed
		if(is.matrix(Gx[[1]])){
			## generate initial values of cline parameters
			initf<-function(L=dim(Gx[[1]])[2],chain_id=1){
        			list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
			dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],
			GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,
			ploidy=pldat,init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_sdk_gl_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
			cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))

		} else{ ## one snp
			## generate initial values of cline parameters
			initf<-function(chain_id=1){
        			list(center=runif(1,.3,.7),v=runif(1,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))		
			dat<-list(L=1,N=length(Gx[[1]]),GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
			H=H,P0=p0,P1=p1,sc=SDc,sv=SDv,ploidy=pldat,init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_one_gl_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	}  else if(hier==TRUE & model=="ancestry" & ploidy=="mixed" & estMu==FALSE){
	## hierachical model, ancestry, mixed ancestry
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,ploidy=pldat,sd0=sd0,init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_z_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
	}  else if(hier==TRUE & model=="ancestry" & ploidy=="mixed" & estMu==TRUE){
	## hierachical model, ancestry, mixed ancestry, estimate means
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,ploidy=pldat,sd0=sd0,mu0=mu0,init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_z_mix_mu,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sv")[[1]],probs=c(.5,.05,.95))
		muc<-quantile(rstan::extract(fit,"muc")[[1]],probs=c(.5,.05,.95))
		muv<-quantile(rstan::extract(fit,"muv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,Muc=muc,Muv=muv,gencline_hmc=fit)
	} else if(hier==FALSE & model=="ancestry" & ploidy=="mixed"){
	## known SD model, ancestry, diploids
		if(is.matrix(Gx)){
			## generate initial values of cline parameters
			initf<-function(L=dim(Gx)[2],chain_id=1){
        			list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
			dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,sc=SDc,sv=SDv,
			ploidy=pldat,init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_sdk_z_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
			cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))

		} else{ ## one snp
			## generate initial values of cline parameters
			initf<-function(chain_id=1){
        			list(center=runif(1,.3,.7),v=runif(1,.9,1.1),alpha=chain_id)
			}
			init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))		
			dat<-list(L=1,N=length(Gx),Z=Gx,H=H,sc=SDc,sv=SDv,ploidy=pldat,
			init=init_ll)
			fit<-rstan::sampling(stanmodels$gencline_one_z_mix,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	}
	
	
	return(Cout)
}
