#' Impose hard sum-to-zero constraints on cline estimates 
#'
#' Re-calculates cline parameter estimates to ensure that the average cline parameter corresponds with expectations for genome-average introgression. 
#' @param hmc HMC object from est_gencline.
#' @param center vector or matrix of cline centers (from est_gencline). If a vector, there should be one element per locus; if a matrix, there should be one row per loucs and one column each for the point estimate and lower and upper bounds of the credible interval.
#' @param v vector or matrix of cline graidents (from est_gencline). If a vector, there should be one element per locus; if a matrix, there should be one row per loucs and one column each for the point estimate and lower and upper bounds of the credible interval.
#' @param transform Boolean variable indicating whether to apply the constraint on the natural scale of v and center (FALSE) or on the transformed scale of log(v) and logit(center)
#' @param ci size of the credible interval, specifically, the equal-tail probability interval, to generate from the HMC object (if supplied) [default = 0.95].
#'
#' @details Genomic clines model locus-specifc patterns of introgression relative to genome-average introgression or ancestry. Thus, if the same loci are used to estimate hybrid index and genomic clines (or if the former are a random sample of the latter), we would expect the average deviations from genome-average admixture to be 0 across loci (i.e., the deviations should cancel out). This is suggested by prior structure of the hierarchical Bayesian model, where the mean for log10(v) and logit(center) are set to 0. This is a soft sum-to-zero constraint (the prior pulls the sum towards zero, but a sum-to-zero constraint is not enforced). This function instead enforced a hard sum-to-zero constraint. This is done either by working with the full HMC output or parameter estimates (just point estiamtes or point estimates and credibel intervals). In the fomrer case, the mean cline parameters (on the log10 or logit scale, as appropriate) at each HMC iteration are forced to be 0 by subtracting off the mean. Point estimates and credible intervals are then re-calculated. In the case of parameter estimates, the mean point estimate (on the log10 or logit scale) is subtracted from the point estimate and credible intervals. Using the full HMC object is generally preferable.
#' @details As an alternative, the constraint can be applied on the natural scale rather than the log or logit scale, that is, the mean cline center can be constrained to 0.5 and the mean gradient (v) to 1. For center, the difference will often be trivial. For gradient the difference could be greater, and my current suggestion is to use the log10 scale as v is a ratio. 
#'
#' @return A list with two vectors (if only a parameter estimate vector was provided) or a matrixes (if a matrix of the full HCM object was given) with the re-calculated, constrained parameter estimates. If the HMC object was given, a point estimate and the bounds of hte specified credible intervals are given.
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC output object.
#'
#' @references
#' Gompert Z, et al. 2024. Bayesian hybrid zone analyses with Hamiltonian Monte Carlo in R. Manuscript in preparation
#' @export
sum2zero<-function (center=NULL, v=NULL, hmc=NULL, transform=TRUE, ci=0.95){ 
	if ((is.null(center) == TRUE | is.null(v) == TRUE) & is.null(hmc)) 
		stop("error, input data are not provided")

	## use the hmc data
	if(is.null(hmc)==FALSE){
		cc<-rstan::extract(fit,"center")[[1]] ## dimensions: rows = HMC steps, cols = loci
		cv<-rstan::extract(fit,"v")[[1]]
		if(transform==TRUE){
			## transform the parameter values and impose sum-to-zero constraints
			t_cv<-log10(cv)
			t_cc<-log(cc/(1-cc))
			for(i in 1:dim(t_cv)[1]){
				t_cv[i,]<-t_cv[i,]-mean(t_cv[i,])
				t_cc[i,]<-t_cc[i,]-mean(t_cc[i,])
			}
			lq<-(1-ci)/2
			uq<-1-lq
			est_v<-apply(t_cv,2,quantile,probs=c(.5,lq,uq))
			est_center<-apply(t_cc,2,quantile,probs=c(.5,lq,uq))
			out<-list(center=est_center,v=est_v)
		} else{
			## don't transform, center should have mean 0.5, gradien mean 1
			for(i in 1:dim(t_cv)[1]){
				cv[i,]<-cv[i,]-mean(cv[i,]) + 1
				cc[i,]<-cc[i,]-mean(cc[i,]) + 0.5
			}
			lq<-(1-ci)/2
			uq<-1-lq
			est_v<-apply(cv,2,quantile,probs=c(.5,lq,uq))
			est_center<-apply(cc,2,quantile,probs=c(.5,lq,uq))
			out<-list(center=est_center,v=est_v)
		}
	} else {
	## use the provide parameter estimates
		# transform parameters and impose sum-to-zero
		if(transform==TRUE){
			t_v<-as.matrix(log10(v))
			t_center<-as.matrix(log(center/(1-center)))
			mns_v<-apply(t_v,2,mean)
			mns_center<-apply(t_center,2,mean)
			for(k in 1:length(mns_v)){
				t_v[,k]<-t_v[,k]-mns_v[1]## subtract off point estimate means
			}
			for(k in 1:length(mns_center)){
				t_center[,k]<-t_center[,k]-mns_center[1]
			}
			out<-list(center=t_center,v=t_v)
		} else{
		## don't transform, center should have mean 0.5, gradien mean 1	
			mns_v<-apply(v,2,mean)
			mns_center<-apply(center,2,mean)
			for(k in 1:length(mns_v)){
				v[,k]<-v[,k]-mns_v[1]+1## subtract off point estimate means
			}
			for(k in 1:length(mns_center)){
				center[,k]<-center[,k]-mns_center[1]+0.5
			}

		}

	}
	return(out)
}

