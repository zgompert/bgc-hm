functions {
	/* computes likelihood weighted by genotype likelihoods */
	real calc_lik(real gl0, real gl1, real gl2, real p, real pl){
		real prob;
		if(pl==2){ /* diploid locus x ind */
            		prob = log(gl0) + log(1-p) + log(1-p); 
            		prob = log_sum_exp(prob, log(gl1) + log(2) + log(p) + log(1-p) + log(2));                  
            		prob = log_sum_exp(prob, log(gl2) + log(p) + log(p));
		} else if(pl==1){ /* haploid locus x ind */
		    prob = log(gl0) + log(1-p);             
            	    prob = log_sum_exp(prob, log(gl1) + log(p));
		} else{
			prob=0;
		}	
		return prob;
	}
}
data{
	int L; /* # of loci */
	int N; /* # of organisms */
	int J; /* # of populations */	
	array[N, L] real<lower=0, upper=1> GL0; /* 2D array of genlik 0*/
	array[N, L] real<lower=0, upper=1> GL1; /* 2D array of genlik 1*/
	array[N, L] real<lower=0, upper=1> GL2; /* 2D array of genlik 2*/
	array[N] real<lower=0, upper=2> ploidy; /* 1D array of ploidy */	
	array[N] int<lower=0> pids; /* 1D array of pop ids for inds. */
}

parameters{
	array[J, L] real<lower=0, upper=1> P; /* 2D array of allele frequencies */
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

