## function to estimate Q, ancestry class model
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
