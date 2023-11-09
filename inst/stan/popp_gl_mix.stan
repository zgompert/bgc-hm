functions {
	/* computes likelihood weighted by genotype likelihoods */
	real calc_lik(real gl0, real gl1, real gl2, real p, real pl){
		real prob;
		if(pl==2){ /* diploid locus x ind */
            		prob = log(gl0) + log(1-p) + log(1-p); 
            		prob = log_sum_exp(prob, log(gl1) + log(2) + log(p) + log(1-p) + log(2));                  
            		prob = log_sum_exp(prob, log(gl2) + log(p) + log(p));
		} else{ /* haploid locus x ind */
		    prob = log(gl0) + log(1-p);             
            prob = log_sum_exp(prob, log(gl2) + log(p));
		}
		return prob;
	}
}
data{
	int L; /* # of loci */
	int N; /* # of organisms */
	int J; /* # of populations */	
	real<lower=0, upper=1> GL0[N, L]; /* 2D array of genlik 0*/
	real<lower=0, upper=1> GL1[N, L]; /* 2D array of genlik 1*/
	real<lower=0, upper=1> GL2[N, L]; /* 2D array of genlik 2*/
	real<lower=0, upper=2> ploidy[N]; /* 1D array of ploidy */	
	int<lower=0> pids[N]; /* 1D array of pop ids for inds. */
}

parameters{
	real<lower=0, upper=1> P[J,L]; /* 2D array of allele frequencies */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(GL0[j,i], GL1[j,i], GL2[j,i], P[pids[j],i], ploidy[j]);
		}
	}
	for(i in 1:L){
	    for(j in 1:J){
    		/* increment priors on allele frequencies */
    		target += beta_lpdf(P[j,i] | 0.5, 0.5);
    	}
    }
}

