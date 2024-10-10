#' Function to calcualte posterior probabilities of extreme cline parameter values
#'
#' Computes the posterior probability that cline parameter values exceed a user defined null or reference value.
#' @param hmc list object produced by est_genocl or est_geocl.
#' @param param character or string specifying the cline parameter to summarize.
#' @param greater Boolean indicating whether to compute the posterior probability the value is (i) greater than (TRUE) or (ii) less than (FALSE) the null or reference value.
#' @param refval null or reference value for posterior probability calculation.
#' 
#' @details This function computes the posterior probability the a parameter value exceeds (or is less than) a null or reference value using samples from the parameter posterior distribution. This works with geographic or genomic cline output and for any model parameter in the HMC output based on the name of the parameter in the HMC object (i.e., Stan model). For genomic clines, this includes cline centers (center), slopes (v), cline standard deviations (sc and sv) and cline means (muc and muv) (not all parameters are components of all models). For geographic clines, this includes  cline centers (cent), slopes (slope), cline widths (w), the mean cline width (mu), and the standard deviation in cline slopes (sigma). In both cases, the name of the parameter in the HMC object is given in parentheses and should be provided as a character or string to the param argument. Different reference or null values are relevant for different questions or parameters. For example, the expected average genomic cline slope, v, is 1, and thus one could identify loci with steeper (or shallower) genomic clines than the genome-average by computing the posterior probability that the slope, v, exceeds (or is less than) 1.
#'
#' @return A vector of posterior probabilities (of length > 1 for vector parameters, such as cline slopes with multiple loci).
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC output object.
#'
#' @references
#' Gompert Z, DeRaad D, Buerkle CA. A next generation of hierarchical Bayesian analyses of hybrid zones enables model-based quantification of variation in introgression in R. bioRxiv 2024.03.29.587395.

#' @export
calc_pp<-function(hmc=NULL,param=NULL,greater=TRUE,refval=NULL){

	## determine whether the hmc object contains gencline_hmc or geo_hmc
	if("gencline_hmc" %in% names(hmc)){
		val<-rstan::extract(hmc$gencline_hmc,param)[[1]]
	} else if("geo_hmc" %in% names(hmc)){
		val<-rstan::extract(hmc$geo_hmc,param)[[1]]
	} else{
		stop("hmc object does not contains output from est_genocl or est_geocl")
	}
	if(is.matrix(val)==TRUE){ ## matrix, multiple parameters
		if(greater==TRUE){
			pp<-apply(val > refval,2,mean)
		} else{
			pp<-apply(val < refval,2,mean)
		}
	} else{ ## just one parameter
		if(greater==TRUE){
			pp<-mean(val > refval)
		} else{
			pp<-mean(val < refval)
		}
	}
	return(pp)
}
