#' Function to estimate geographic clines for a set of genetic loci or hybrid index
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of geographic clines from genetic data. 
#' @param G genetic data for the sampled hybrid zone in the form of a matrix for known genotypes (rows = individuals, columns = loci), a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype), or  matrix of ancestry for the ancestry model (rows = individuals, columns = loci).
#' @param P matrix of allele frequencies for the hybrid zone with rows denoting sampled populations (localities) and columns denoting loci.
#' @param Geo vector of geographic coordinates for the sampled populations. 
#' @param Ids indexes designating which population each individual belongs to. This should be provided as a single vector (length equal to the number of sampled individuals). The indexes should range from 1 to the number of populations.
#' @param model for genetic data, either 'genotype' for known genotypes or 'glik' for genotype likelihoods.
#' @param ploidy specifies ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.
#' @param pldat matrix or list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid).
#' @param hier Boolean, fit hierarchical model (TRUE) that estimates cline slope SDs or non-hierarchical model for one locus with (relatively) uninformative priors (FALSE).
#' @param prec approximate precision for known (or estimated) allele frequencies, which should be set to about 1/2N (use the average 2N across populations). Do not set this to 0; this can result in errors when working with the log of the allele frequencies.
#' @param y_lb minimum value of logit allele frequencies to include in the model (the default of -2 is a good choice).
#' @param y_ub minimum value of logit allele frequencies to include in the model (the default of 2 is a good choice).
#' @param gamma_a alpha parameter for gamma priors on the error standard deviation and standard deviation on the prior for slope (default = 0.1).
#' @param gamma_b beta parameter for gamma priors on the error standard deviation and standard deviation on the prior for slope (default = 0.01).
#' @param n_chains number of HMC chains for posterior inference.
#' @param n_iters A positive integer specifying the number of iterations for each chain (including warmup), default is 2000.
#' @param p_warmup proportion (between 0 and 1) of n_iters to use as warmup (i.e., burnin), default is 0.5.
#' @param n_thin positive integer, save every n_thin HMC iterations for the posterior, default is 1.
#' @param n_cores number of cores to use for HMC, leave as NULL to automatically detect the number of cores present (no more than n_chains cores can be used even if more cores are available).
#' @param full boolean denoting whether (TRUE, the default) or not (FALSE) to return the HMC Stan object with the full set of samples from the posterior.
#'
#' @details
#' Geographic clines are estimated from population (deme) allele frequencies. This is done using a linear model for the logit of the allele frequencies (a sigmoid cline on the natural scale becomes linear on the logit scale). Users can provide allele frequency estimates or the allele frequencies can be estimate from genotypic data. In the latter case, allele frequencies are first estimated based on known genotypes (model = 'genotype') or genotype likelihoods (model = 'glik') using an analytical solution for the posterior. Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, use 0 and 1. Genotype likelihoods should be on their natural scale (not phred scaled) and the values for each locus for an individual should sum to 1 (i.e., the likelihoods are scaled to be probabilities). The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, use the 0 and 1 matrixes to store the likelihoods of the two possible states. If provided directly, allele frequencies should be given as a matrix, with one column per locus (assumes bi-allelic SNPs or the equivalent) and one row per population (deme).  
#' @details
#' Ploidy data are only required for mixed ploidy genotypic data (they are not used if allele frequencies are provided directly). In this case, there should be one matrix with the same dimensions as the genetic data. The values in the matrix indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1). Ploidy can be set to 0 to denote missing data.
#' @details
#' The model assumes organisms have been sampled from populations (demes) along a 1D transect through a hybrid zone. Various approaches exist for approximating a 2D sampling scheme in 1D and can be used to transform coordinates to a single dimension. No specific coordinate units are expected, and coordinates are always centered (given a mean of 0) prior to analysis (this is done internally if not done by the user). If population allele frequencies are given directly, populations are assumed to be in the same order in the Geo vector and allele frequency matrix. If genotypic data are provided, and additional object, Ids, is required that indicates which population (numbered 1 to the number of populations and following the order in Geo) each individual belongs to. 
#' @details The model works with the logit of the allele frequencies. Consequently, allele frequencies of 0 are not allowed (these will cause an error). The value specified by prec will be added to allele frequencies of 0 and subtracted from allele frequencies of 1. This prevents problems with taking logs and also is meant to reflect that fact that one cannot be certain an allele is not present in a population. We recommend setting prec to 1/2N, where N is the mean (or median) sample size across demes. Single and multilocus clines should be well approximated by a linear function for the logit allele frequencies near the center of the hybrid zone; this is the reason for only analyzing populations with intermediate allele frequencies (logit p between y_lb and y_ub, -2 and 2 by default, or p of about 0.11 to 0.88)
#'
#' @return A list of parameter estimates and full HMC results from stan. Estimates are provided for the cline width on the natural scale for each locus (w), the slope on the logit scale (slope), the cline center based on the centered geographic coordinates on the logit scale (center), and the mean (mu) and standard deviation (sigma) for cline widths on the logit scale (the latter two only apply to the hierarchical model). Parameter estimates are provided as a point estimate (median of the posterior), 90% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution), and 95% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. The full HMC output from rstan is provided as the final element in the list. This can be used for HMC diagnostics and to extract other model outputs not provided by default. 
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC output object.
#'
#' @references
#' Gompert Z, DeRaad D, Buerkle CA. A next generation of hierarchical Bayesian analyses of hybrid zones enables model-based quantification of variation in introgression in R. bioRxiv 2024.03.29.587395.

#' @export
#' @examples
#' \dontrun{
#' ## load the data set
#' data(pfreqs)
#' ## this includes one object, a matrix P with allele frequencies
#' ## 110 rows = demes, 51 columns = loci
#'
#' ## use standardized deme numbers as geographic coordinates
#' x<-1:110
#' geo<-(x-mean(x))/sd(x)
#'
#' ## fit the geographic cline model
#' o<-est_geocl(P=P,Geo=geo,prec=0.01,y_lb=-2,y_ub=2,hier=TRUE,n_iters=5000)
#'
#' ## plot clines on logit scale, which should be linear
#' plot(geo,o$cent[1,1] + o$slope[1,1] * geo,type='l',ylim=c(-15,15),ylab="Logit allele frequency",xlab="Deme number",
#'	axes=FALSE)
#' axis(1,at=geo[seq(5,110,5)],x[seq(5,110,5)])
#' axis(2)
#' box()
#' for(i in 2:51){
#'	lines(geo,o$cent[i,1] + o$slope[i,1] * geo)
#' }
#' }
est_geocl<-function(G=NULL,P=NULL,Geo=NULL,Ids=NULL,model="genotype",ploidy="diploid",pldat=NULL,hier=TRUE,prec=0.001,
	y_lb=-2,y_ub=2,gamma_a=0.1,gamma_b=0.01,n_chains=4,n_iters=2000,p_warmup=0.5,n_thin=1,n_cores=NULL,
	full=TRUE){

        ## get or set number of cores for HMC
        if(is.null(n_cores)){
                options(mc.cores = parallel::detectCores())
        } else{
                mc.cores<-n_cores
        }
        
        ## determine number of warmup iterations
        n_warmup<-floor(p_warmup * n_iters)
	
	if(is.null(P)){ 
     	## need to first estimate population allele frequencies
       	## use analytical solution for the posterior
		if(model=="genotype" & ploidy=="diploid"){
			P<-matrix(NA,nrow=length(Geo),ncol=dim(G)[2])
			for(i in 1:dim(G)[2]){
				Y<-tapply(X=G[,i],INDEX=Ids,sum)
				N<-tapply(X=G[,i]>=0,INDEX=Ids,sum)
				# posterior median
				P[,i]<-qbeta(0.5,Y+0.5,N-Y+0.5)
			}
		} else if(model=="genotype" & ploidy=="mixed"){
			P<-matrix(NA,nrow=length(Geo),ncol=dim(G)[2])
			for(i in 1:dim(G)[2]){
				Y<-tapply(X=G[,i],INDEX=Ids,sum)
				N<-tapply(X=pldat[,i],INDEX=Ids,sum)
				# posterior median
				P[,i]<-qbeta(0.5,Y+0.5,N-Y+0.5)
			}
		} else if(model=="glik" & ploidy=="diploid"){
			P<-matrix(NA,nrow=length(Geo),ncol=dim(G)[2])
			for(i in 1:dim(G)[2]){
				Y<-tapply(X=G[[2]][,i]+G[[3]][,i]*2,INDEX=Ids,sum)
				N<-tapply(X=G[[1]][,i]>=0,INDEX=Ids,sum)
				# posterior median
				P[,i]<-qbeta(0.5,Y+0.5,N-Y+0.5)
			}		
		} else { ## glik and mixed
			P<-matrix(NA,nrow=length(Geo),ncol=dim(G)[2])
			for(i in 1:dim(G)[2]){
				Y<-tapply(X=G[[2]][,i]+G[[3]][,i]*2,INDEX=Ids,sum)
				N<-tapply(X=pldat[,i],INDEX=Ids,sum)
				# posterior median
				P[,i]<-qbeta(0.5,Y+0.5,N-Y+0.5)
			}
		}
	}
	## avoids log(0) and divisions by 0
	P[P < prec]<-prec
	P[P > (1-prec)]<-1-prec
    	
	## center geography
	Geo<-Geo-mean(Geo)
	Y<-log(P/(1-P))
	if(hier==TRUE){
		dat<-list(L=dim(P)[2],J=dim(P)[1],Y=as.matrix(log(P/(1-P))),geo=Geo,lb=y_lb,ub=y_ub,
			ga=gamma_a,gb=gamma_b)
		fit<-rstan::sampling(stanmodels$geocline,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		w<-t(apply(rstan::extract(fit,"w")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		cent<-t(apply(rstan::extract(fit,"cent")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		slope<-t(apply(rstan::extract(fit,"slope")[[1]],2,quantile,probs=c(.5,.025,.05,.95,.975)))
		mu<-quantile(rstan::extract(fit,"mu")[[1]],probs=c(.5,.025,.05,.95,.975))
		sigma<-quantile(rstan::extract(fit,"sigma")[[1]],probs=c(.5,.025,.05,.95,.975))
		## create a list with parameter estimates plus full hmc object
		if(full==TRUE){
			geoout<-list(w=w,cent=cent,slope=slope,mu=mu,sigma=sigma,geo_hmc=fit)
		} else{
			geoout<-list(w=w,cent=cent,slope=slope,mu=mu,sigma=sigma)
		}	
	} else {
		P<-as.vector(P)
		dat<-list(J=length(P),Y=(log(P/(1-P))),geo=Geo,lb=y_lb,ub=y_ub,
			ga=gamma_a,gb=gamma_b)
		fit<-rstan::sampling(stanmodels$geocline_one,data=dat,
			iter=n_iters,warmup=n_warmup,thin=n_thin)
		w<-quantile(rstan::extract(fit,"w")[[1]],2,probs=c(.5,.025,.05,.95,.975))
		cent<-quantile(rstan::extract(fit,"cent")[[1]],2,probs=c(.5,.025,.05,.95,.975))
		slope<-quantile(rstan::extract(fit,"slope")[[1]],2,probs=c(.5,.025,.05,.95,.975))
		## create a list with parameter estimates plus full hmc object
		if(full==TRUE){
			geoout<-list(w=w,cent=cent,slope=slope,geo_hmc=fit)
		} else{
			geoout<-list(w=w,cent=cent,slope=slope)
		}
	}
	return(geoout)
}
