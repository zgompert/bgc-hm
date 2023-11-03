#' Estimate genomic clines using a hierarchical model with unknown cline variances 
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of genomic clines from genetic data. This fits a hierarchical log-logistic model for genomic clines with two key parameters, cline center and cline gradient (i.e., slope, inversely proportional to cline width). The model also estimates the cline SDs across loci, SDc (variation in logit centers) and SDv (variation in log10 gradients). 
#' @param Gx genetic data for putative hybrids, matrix for known genotypes, list of matrixes for genotype likelihoods, matrix of ancestry for the ancestry model
#' @param G0 genetic data for parental reference set 0
#' @param G1 genetic data for parental reference set 1
#' @param p0 allele frequencies for parental reference set 0
#' @param p1 allele frequencies for parental reference set 1
#' @param H hybrid indexes for the putative hybrids
#' @param model genotype, glik, or ancestry
#' @param ploidy diploid or mixed
#' @param pldat ploidy data for mixed ploidy
#' @param hier boolean, fit hierarchical model (TRUE) or assume cline SDs known (FALSE)
#' @param SDc known cline center SD on logit scale
#' @param SDv known cline gradient SD on log10 scale
#' @param n_chains number of HMC chains for posterior inference
#'
#' @return list of parameter estimates and full HMC results from stan, this includes cline parameters (center and gradient) and standard deviations describing variability in clines across loci (SDc and SDv) 
#'
#' @export
est_genocl<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,H=NULL,model="genotype",ploidy="diploid",pldat=NULL,
	hier=TRUE,SDc=NULL,SDv=NULL,n_chains=4){
	
	options(mc.cores = parallel::detectCores())
	
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
	}
	
	
	## fit appropriate genomic cline model based on options
	if(hier==TRUE & model=="genotype" & ploidy=="diploid"){
	## hierachical model, known genotypes, diploids
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,
		init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline,data=dat)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
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
			fit<-rstan::sampling(stanmodels$gencline_sdk,data=dat)
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
			fit<-rstan::sampling(stanmodels$gencline_one,data=dat)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)

	} else if(hier==TRUE & model=="glik" & ploidy=="diploid"){
	## hierachical model, genotype likelihoods, diploids
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
		H=H,P0=p0,P1=p1,init=init_ll)

		if(dim(Gx[[1]])[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_gl,data=dat)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
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
			fit<-rstan::sampling(stanmodels$gencline_sdk_gl,data=dat)
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
			fit<-rstan::sampling(stanmodels$gencline_one_gl,data=dat)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="ancestry" & ploidy=="diploid"){
	## hierachical model, ancestry, diploids
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_z,data=dat)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
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
			fit<-rstan::sampling(stanmodels$gencline_sdk_z,data=dat)
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
			fit<-rstan::sampling(stanmodels$gencline_one_z,data=dat)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)

	} else if(hier==TRUE & model=="genotype" & ploidy=="mixed"){
	## hierachical model, known genotypes, mixed ploidy
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,ploidy=pldat,
		init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_mix,data=dat)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
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
			fit<-rstan::sampling(stanmodels$gencline_sdk_mix,data=dat)
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
			fit<-rstan::sampling(stanmodels$gencline_one_mix,data=dat)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	} else if(hier==TRUE & model=="glik" & ploidy=="mixed"){
	## hierachical model, genotype likelihoods, mixed ploidy
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
		H=H,P0=p0,P1=p1,ploidy=pldat,init=init_ll)

		if(dim(Gx[[1]])[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_gl_mix,data=dat)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
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
			fit<-rstan::sampling(stanmodels$gencline_sdk_gl_mix,data=dat)
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
			fit<-rstan::sampling(stanmodels$gencline_one_gl_mix,data=dat)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	}  else if(hier==TRUE & model=="ancestry" & ploidy=="mixed"){
	## hierachical model, ancestry, mixed ancestry
		## generate initial values of cline parameters
		initf<-function(L=dim(Gx)[2],chain_id=1){
        		list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
		}
		init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,H=H,ploidy=pldat,init=init_ll)

		if(dim(Gx)[2]==1){stop("at least two loci required for the hierarchical model")}

		fit<-rstan::sampling(stanmodels$gencline_z_mix,data=dat)
		cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
		cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
		sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
		sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)
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
			fit<-rstan::sampling(stanmodels$gencline_sdk_z_mix,data=dat)
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
			fit<-rstan::sampling(stanmodels$gencline_one_z_mix,data=dat)
			cc<-quantile(rstan::extract(fit,"center")[[1]],probs=c(.5,.05,.95))
			cv<-quantile(rstan::extract(fit,"v")[[1]],probs=c(.5,.05,.95))
		}
		## create a list with parameter estimates plus full hmc object
        	Cout<-list(center=cc,gradient=cv,gencline_hmc=fit)
	}
	
	
	return(Cout)
}
