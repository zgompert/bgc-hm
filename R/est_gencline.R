#' Estimate genomic clines using a hierarchical model with unknown cline variances 
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of genomic clines from genetic data. This fits a hierarchical log-logistic model for genomic clines with two key parameters, cline center and cline gradient (i.e., slope, inversely proportional to cline width). The model also estimates the cline SDs across loci, SDc (variation in logit centers) and SDv (variation in log10 gradients). 
#' @param Gx genetic data for putative hybrids
#' @param G0 genetic data for parental reference set 0
#' @param G1 genetic data for parental reference set 1
#' @param p0 allele frequencies for parental reference set 0
#' @param p1 allele frequencies for parental reference set 1
#' @param H hybrid indexes for the putative hybrids
#' @param model genotype, glik, or ancestry
#' @param ploidy diploid or mixed
#' @param pldat ploidy data for mixed ploidy
#'
#' @return list of parameter estimates and full HMC results from stan, this includes cline parameters (center and gradient) and standard deviations describing variability in clines across loci (SDc and SDv) 
#'
#' @export
est_genocl<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,H=NULL,model="genotype",ploidy="diploid",pldat=NULL){
	
	n_chains<-4 # make this an argument?
	initf<-function(L=dim(Gx)[2],chain_id=1){
        	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
	}
	init_ll<-lapply(1:n_chains, function(id) initf(chain_id = id))

	dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,H=H,P0=p0,P1=p1,init=init_ll)

	fit<-rstan::sampling(stanmodels$gencline,data=dat)
	cc<-t(apply(rstan::extract(fit,"center")[[1]],2,quantile,probs=c(.5,.05,.95)))
	cv<-t(apply(rstan::extract(fit,"v")[[1]],2,quantile,probs=c(.5,.05,.95)))
	sdc<-quantile(rstan::extract(fit,"sdc")[[1]],probs=c(.5,.05,.95))
	sdv<-quantile(rstan::extract(fit,"sdv")[[1]],probs=c(.5,.05,.95))
	## create a list with parameter estimates plus full hmc object
        Cout<-list(center=cc,gradient=cv,SDc=sdc,SDv=sdv,gencline_hmc=fit)

	return(Cout)
}
