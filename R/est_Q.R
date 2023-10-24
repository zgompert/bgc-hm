#' Function to estimate ancestry class proportions (Q)
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of ancestry classes from genetic data. Ancestry classes denote the proportion of an individual's genome where both gene copies come from source 1 (Q11), both gene copies come from source 0 (Q00), or where one gene copy comes from source 1 and one from source 0.
#' @param Gx genetic data for putative hybrids
#' @param G0 genetic data for parental reference set 0
#' @param G1 genetic data for parental reference set 1
#' @param p0 allele frequencies for parental reference set 0
#' @param p1 allele frequencies for parental reference set 1
#' @param model genotype, glik, or ancestry
#' @param ploidy diploid or mixed
#' @param pldat ploidy data for mixed ploidy
#'
#' @return list of parameter estimates and full HMC results from stan, this includes Q (ancestry class proportions) and hybrid index, which is derived from Q
#'
#' @export
est_Q<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,model="genotype",ploidy="diploid",pldat=NULL){

	dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1)
	fit<-rstan::sampling(stanmodels$Q,data=dat)
	Q<-extract(fit,"Q")[[1]]
	Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
	Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
	Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
	hi<-t(apply(extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
	## create a list with parameter estimates plus full hmc object
        Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi,Q_hmc=fit)

	return(Qout)
}
