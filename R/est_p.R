#' Function to estimate parental allele frequencies
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference parental allele frequencies from genetic data.
#' @param Gx genetic data for putative hybrids
#' @param G0 genetic data for parental reference set 0
#' @param G1 genetic data for parental reference set 1
#' @param model genotype, glik, or ancestry
#' @param ploidy diploid or mixed
#' @param pldat ploidy data for mixed ploidy
#'
#' @return list of parameter estimates and full HMC results from stan
#'
#' @export
est_p<-function(G0=NULL,G1=NULL,model="genotype",ploidy="diploid",pldat=NULL){

	options(mc.cores = parallel::detectCores())

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
