#' Function to estimate ancestry class proportions (Q) using a variational algorithm for approximate posterior sampling
#'
#' Uses Stan's variational algorithm for approximate posterior sampling for Bayesian inference  of ancestry classes proportions from genetic data. Ancestry classes denote the proportion of an individual's genome where both gene copies come from source 1 (Q11), both gene copies come from source 0 (Q00), or where one gene copy comes from source 1 and one from source 0 (Q10).This is an approximation to the full Bayesian analysis with HMC. This is an experimental feature and might not perform well in all cases, but will generally be much faster than HMC.
#' @param Gx genetic data for putative hybrids in the form of a matrix for known genotypes (rows = individuals, columns = loci), a list of matrixes for genotype likelihoods (same dimensions but one matrix per genotype), or  matrix of ancestry for the ancestry model (rows = individuals, columns = loci).
#' @param G0 genetic data for parental reference set 0 formatted as described for Gx.
#' @param G1 genetic data for parental reference set 1 formatted as described for Gx.
#' @param p0 vector of allele frequencies for parental reference set 0 (one entry per locus).
#' @param p1 vector allele frequencies for parental reference set 1 (one entry per locus).
#' @param model for genetic data, either 'genotype' for known genotypes, 'glik' for genotype likelihoods, or 'ancestry' for known ancestry.
#' @param ploidy species ploidy, either all 'diploid' or 'mixed' for diploid and haploid loci or individuals.
#' @param pldat matrix or list of matrixes of ploidy data for mixed ploidy (rows = individuals, columns = loci) indicating ploidy (2 = diploid, 1 = haploid).
#' @param method variational inference algorithm to use, either 'meanfield’ or 'fullrank’, default is meanfiled; fullrank deals better with covariances among parameters in the (approximate and transformed) posterior distribution
#'
#' @return A list of parameter estimates, this includes Q (ancestry class proportions) and hybrid indexes, which are derived from Q. Parameter estimates are provided as a point estimate (median of the posterior) and 95% equal-tail probability intervals (2.5th and 97.5th quantiles of the posterior distribution). These are provided as a vector or matrix depending on the dimensionality of the parameter. 
#'
#' @details Ancestry class proportions can be estimated from known genotypes (model = 'genotype'), genotype likelihoods, (model = 'glik') or known (estimated) ancestry (model = 'ancestry'). Genotypes should be encoded as 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote). No specific polarization (e.g., minor allele, reference allele, etc.) of 0 vs 2 is required. For haploid loci, use 0 and 1. Genotype likelihoods should be on their natural scale (not phred scaled) and the values for each locus for an individual should sum to 1. The data should be provided as a list of three matrixes, with the matrixes giving the likelihoods for genotypes 0, 1 and 2 respectively. Thus, each matrix will have one row per individual and one column per locus. For haploid loci with genotype likelihoods, you must use the 0 and 1 matrixes to store the likelihoods of the two possible states. For the ancestry model, hybrid indexes are inferred directly from known local (locus-specific) ancestry rather than from genotype data. Users are free to use whatever software they prefer for local ancestry inference (many exist). In this case, each entry in the individual (rows) by locus (columns) matrix should denote the number of gene copies inherited from parental population 1 (where pure parent 1 corresponds with a hybrid index of 1 and pure parent 0 corresponds with a hybrid index of 0). Haploids can be encoded using 0 and 1.  For all models, missing data can be encoded by setting the ploidy for an individual/locus to 0 (this indicates no information, whereas genotype likelihoods encode uncertainty in genotypes) and the genotype to NA. 
#' @details
#' Hybrid genetic (or ancestry) data are always required. For genotype or genotype likelihood models, users must either provide pre-estimated parental allele frequencies or parent genetic (genotypes or genotype likelihoods) that can be used to infer allele frequencies. Parental data are not required for the ancestry model.
#' @details
#' Ploidy data are only required for the mixed ploidy data. In this case, there should be one matrix for the hybrids or a list of matrixes for the hybrids (1st matrix) and each parent (2nd and 3rd matrixes, with parent 0 first). The latter is required for the genotype or genotype likelihood models if parental allele frequencies are not provided. The matrixes indicate whether each locus (column) for each individual (row) is diploid (2) or haploid (1) (use 0 for missing data). Haploid loci are useful for estimating the proportion of the genome inherited from each reference population and thus are relevant for this model, but do not specifically provide information about how this is partitition into heterozygous vs homozygous ancestry classes.
#'
#' @seealso 'rstan::stan' for details on HMC with stan and the rstan HMC output object.
#'
#' @references
#' Gompert Z, DeRaad D, Buerkle CA. 2024. A next generation of hierarchical Bayesian analyses of hybrid zones enables model-based quantification of variation in introgression in R. Ecology and Evolution, 14:e70584.
#'
#' Kucukelbir A, Tran D, Ranganath R, Gelman A, Blei DM. 2017. Automatic differentiation variational inference. Journal of Machine Learning Research, 18:1-45.
#' @export
#' @examples
#'\dontrun{
#' ## load the data set
#' data(genotypes)
#' ## this includes three objects, GenHybrids, GenP0, and GenP1

#' ## estimate parental allele frequencies, uses default HMC settings
#' p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")
#' 
#' ## estimate interspecific ancestry, this can
#' ## be especially informative about the types of hybrids present
#' q_out<-est_Q_vi(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")
#'}

est_Q_vi<-function(Gx=NULL,G0=NULL,G1=NULL,p0=NULL,p1=NULL,model="genotype",ploidy="diploid",pldat=NULL,
		method="meanfield"){

        if(is.list(Gx)==TRUE & model!="glik"){
                stop("List input is only valide for the glik model")
        } else if(is.list(Gx)==FALSE & model=="glik"){
                stop("List input required for the genotype likelihood model")
        } else if(! model %in% c("genotype","glik","ancestry")){
            stop("Unknown model, must be genotype, glik, or ancestry")
        } else if(! ploidy %in% c("diploid","mixed")){
            stop("Only diploid and mixed are accepted for the ploidy argument")
        } else if(ploidy=="mixed" & is.null(pldat)){
            stop("pldat cannot be NULL for mixed ploitdy")
        }

        ## estimate parental allele frequencies if not provided
        if(is.null(p0) | is.null(p1)){
                po<-est_p(G0=G0,G1=G1,model=model,ploidy=ploidy,pldat=pldat,HMC=FALSE)
                p0<-po$p0[,1]
                p1<-po$p1[,1]
        }

	if(sum(p0==0&p1==0)+sum(p0==1&p1==1) > 0){
		message("one or more loci are fixed for the same allele in both parental source populations, this can cause problems\n")
	}

        ## now have population allele frequencies, just need hybrid ploidy data if it
        ## was provided as a list
        if(is.list(pldat)==TRUE){
                pldat<-pldat[[1]]
        }

        ## fit the appropriate hybrid index model based on options
        if(model=="genotype" & ploidy=="diploid"){
                ## diploid with known genotypes
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1)
		fit<-rstan::vb(stanmodels$Q,data=dat,algorithm=method)
		Q<-rstan::extract(fit,"Q")[[1]]
		Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
		Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
		Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
		## create a list with parameter estimates
		Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi)
	} else if(model=="glik" & ploidy=="diploid"){
                ## diploid with genotype likelihoods
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]]
			  ,P0=p0,P1=p1)
		fit<-rstan::vb(stanmodels$Q_gl,data=dat,algorithm=method)
		Q<-rstan::extract(fit,"Q")[[1]]
		Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
		Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
		Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
		## create a list with parameter estimates
		Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi)
	} else if(model=="ancestry" & ploidy=="diploid"){
                ## diploid with known ancestry
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx)
		fit<-rstan::vb(stanmodels$Q_z,data=dat,algorithm=method)
		Q<-rstan::extract(fit,"Q")[[1]]
		Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
		Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
		Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
		## create a list with parameter estimates
		Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi)
	} else if(model=="genotype" & ploidy=="mixed"){
        ## mixed ploidy with known genotypes
        Gx[pldat==0]<-0
		if(sum(is.na(Gx)) > 0){ 
		    stop("All NA genotypes must have ploidy set to 0")
		}            
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],G=Gx,P0=p0,P1=p1,ploidy=pldat)
		fit<-rstan::vb(stanmodels$Q_mix,data=dat,algorithm=method)
		Q<-rstan::extract(fit,"Q")[[1]]
		Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
		Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
		Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
		## create a list with parameter estimates
		Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi)
	} else if(model=="glik" & ploidy=="mixed"){
                ## mixed ploidy with known genotypes
		dat<-list(L=dim(Gx[[1]])[2],N=dim(Gx[[1]])[1],GL0=Gx[[1]],GL1=Gx[[2]],GL2=Gx[[3]],
			  P0=p0,P1=p1,pl=pldat)
		fit<-rstan::vb(stanmodels$Q_gl_mix,data=dat,algorithm=method)
		Q<-rstan::extract(fit,"Q")[[1]]
		Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
		Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
		Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
		## create a list with parameter estimates
		Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi)
	} else if(model=="ancestry" & ploidy=="mixed"){
                ## mixed ploidy with known ancestry
		dat<-list(L=dim(Gx)[2],N=dim(Gx)[1],Z=Gx,pl=pldat)
		fit<-rstan::vb(stanmodels$Q_z_mix,data=dat,algorithm=method)
		Q<-rstan::extract(fit,"Q")[[1]]
		Q11<-t(apply(Q[,,1],2,quantile,probs=c(.5,.05,.95)))
		Q10<-t(apply(Q[,,2],2,quantile,probs=c(.5,.05,.95)))
		Q00<-t(apply(Q[,,3],2,quantile,probs=c(.5,.05,.95)))
		hi<-t(apply(rstan::extract(fit,"H")[[1]],2,quantile,probs=c(.5,.05,.95)))
		## create a list with parameter estimates
		Qout<-list(Q11=Q11,Q10=Q10,Q00=Q00,hi=hi)
	}
	return(Qout)
}
