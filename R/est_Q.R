## function to estimate hybrid index
est_hi<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,model="genotype",ploidy="diploid",pldat=NULL){

	dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1)
	fit<-rstan::sampling(stanmodels$hi,data=dat)
	hi<-t(apply(extract(fit,"H")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
	## create a list with parameter estimates plus full hmc object
        hiout<-list(hi=hi,hi_hmc=fit)

	return(hiout)
}
