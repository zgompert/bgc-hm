#' Function to estimate hybrid index
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of hybrid index from genetic data.
#' @param Gx genetic data for putative hybrids
#' @param G0 genetic data for parental reference set 0
#' @param G1 genetic data for parental reference set 1
#' @param p0 allele frequencies for parental reference set 0
#' @param p1 allele frequencies for parental reference set 1
#' @param model genotype, glik, or ancestry
#' @param ploidy diploid or mixed
#' @param pldat ploidy data for mixed ploidy
#'
#' @return list of parameter estimates and full HMC results from stan
#'
#' @export
est_hi<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,model="genotype",ploidy="diploid",pldat=NULL){

	dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1)
	fit<-rstan::sampling(stanmodels$hi,data=dat)
	hi<-t(apply(extract(fit,"H")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
	## create a list with parameter estimates plus full hmc object
        hiout<-list(hi=hi,hi_hmc=fit)

	return(hiout)
}
