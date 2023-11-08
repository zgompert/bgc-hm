#' Function to estimate geographic clines for a set of genetic loci or hybrid index
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of geographic clines from genetic data. 
#' @param G genetic data 
#' @param P allele frequency data
#' @param Geo geographic coordinates for the sampled populations, these will be centered automatically
#' @param Ids indexes designating which population each individual belongs to, indexes run from 1 to the number of populations
#' @param model genotype or glik
#' @param ploidy diploid or mixed
#' @param pldat ploidy data for mixed ploidy
#' @param prec precision for known allele frequencies, set to approximately 1/2N
#'
#' @return list of parameter estimates and full HMC results from stan
#'
#' @export
est_geocl<-function(G=NULL,P=NULL,Geo=NULL,Ids=NULL,model="genotype",ploidy="diploid",pldat=NULL,prec=0.001){

	options(mc.cores = parallel::detectCores())
        
    if(is.null(P)){ 
     	## need to first estimate population allele frequencies
       	## use analytical solution for the posterior
       	if(model=="genotype" & ploidy=="diploid"){
       	    P<-matrix(NA,nrow=length(Geo),ncol=dim(G)[2])
       	    for(i in 1:dim(G)[2]){
       		    Y<-tapply(X=G[,i],INDEX=Ids,sum)
       		    N<-tapply(X=G[,i]>=0,INDEX=Ids,sum)
       		    P[,i]<-pbeta(q=0.5,Y+0.5,N-Y+0.5)
       	    }
	    } else if(model=="genotype" & ploidy=="mixed"){
	        P<-matrix(NA,nrow=length(Geo),ncol=dim(G)[2])
       	    for(i in 1:dim(G)[2]){
       		    Y<-tapply(X=G[,i],INDEX=Ids,sum)
       		    N<-tapply(X=pldat[,i],INDEX=Ids,sum)
       		    P[,i]<-pbeta(q=0.5,Y+0.5,N-Y+0.5)
       	    }
	    } else if(model=="glik" & ploidy=="diploid"){
	        dat<-list(L=dim(G[[1]])[2],N=dim(G[[1]])[1],J=length(Geo),GL0=G[[1]],GL1=G[[2]],GL2=G[[3]],pids=Ids)
	        fit<-rstan::sampling(stanmodels$popp_gl,data=dat)
	        Px<-rstan::extract(fit,"P")[[1]]
	        P<-apply(Px,c(2,3),median) 
	    } else { ## glik and mixed
	        dat<-list(L=dim(G[[1]])[2],N=dim(G[[1]])[1],J=length(Geo),GL0=G[[1]],GL1=G[[2]],GL2=G[[3]],pids=Ids)
	        fit<-rstan::sampling(stanmodels$popp_gl,data=dat)
	        Px<-rstan::extract(fit,"P")[[1]]
	        P<-apply(Px,c(2,3),median) 
	    }
	}
	## avoids log(0) and divisions by 0
    P[P < prec]<-prec
    P[P > (1-prec)]<-1-prec
    	
	## center geography
	Geo<-Geo-mean(Geo)
		        
	dat<-list(L=dim(P)[2],J=dim(P)[1],P=log(P/(1-P)),geo=Geo)
	fit<-rstan::sampling(stanmodels$geocline,data=dat)
	w<-t(apply(rstan::extract(fit,"w")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
	mu<-quantile(rstan::extract(fit,"mu")[[1]],probs=c(.5,.025,.05,.95,.975))
	sigma<-quantile(rstan::extract(fit,"sigma")[[1]],probs=c(.5,.025,.05,.95,.975))
	## create a list with parameter estimates plus full hmc object
	geoout<-list(w=w,mu=mu,sigma=sigma,geo_hmc=fit)
	return(geoout)
}
