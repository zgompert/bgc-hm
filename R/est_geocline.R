#' Function to estimate geographic clines for a set of genetic loci or hybrid index
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of geographic clines from genetic data. 
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
est_geocl<-function(G=NULL,Geo=NULL,Ids=NULL,model="genotype",ploidy="diploid",pldat=NULL){

   options(mc.cores = parallel::detectCores())
        
	dat<-list(L=dim(G)[2],N=dim(G)[1],J=length(Geo),G=G,geo=Geo,popid=Ids)
	fit<-rstan::sampling(stanmodels$geocline,data=dat)
	w<-t(apply(rstan::extract(fit,"w")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
	mu<-quantile(rstan::extract(fit,"mu")[[1]],probs=c(.5,.025,.05,.95,.975))
    sigma<-quantile(rstan::extract(fit,"sigma")[[1]],probs=c(.5,.025,.05,.95,.975))
	## create a list with parameter estimates plus full hmc object
    geoout<-list(w=w,mu=mu,sigma=sigma,geo_hmc=fit)

	return(geoout)
}
