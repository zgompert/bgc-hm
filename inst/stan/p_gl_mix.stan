functions {
	/* computes likelihood weighted by genotype likelihoods */
	real calc_lik(real gl0, real gl1, real gl2, real p, int pl){
		real prob;
		if(pl==2){ /* diploid locus x ind */
        		prob = log(gl0) + log(1-p) + log(1-p); 
       			prob = log_sum_exp(prob, log(gl1) + log(2) + log(p) + log(1-p));                  
        		prob = log_sum_exp(prob, log(gl2) + log(p) + log(p));
		} else if(pl==1) { /* haploid locus x ind */
			prob = log(gl0) + log(1-p); 
			/* store one copy alt allele as gl1 */
        		prob = log_sum_exp(prob, log(gl1) + log(p));
        	} else {
        		prob = 0; 
        	}	
		return prob;
	}
}
data{
	int L; /* # of loci */
	int N0; /* # of parent 0 organisms */
	int N1; /* # of parent 1 organisms */	
	array[N0, L] real<lower=0, upper=2> GL00; /* 2D array of genlik 0 for parent 0 */
	array[N1, L] real<lower=0, upper=2> GL10; /* 2D array of genlik 0 for parent 1 */
	array[N0, L] real<lower=0, upper=2> GL01; /* 2D array of genlik 1 for parent 0 */
	array[N1, L] real<lower=0, upper=2> GL11; /* 2D array of genlik 1 for parent 1 */
	array[N0, L] real<lower=0, upper=2> GL02; /* 2D array of genlik 2 for parent 0 */
	array[N1, L] real<lower=0, upper=2> GL12; /* 2D array of genlik 2 for parent 1 */
	array[N0, L] int<lower=0, upper=2> ploidy0; /* 2D array of ploidy for parent 0 */
	array[N1, L] int<lower=0, upper=2> ploidy1; /* 2D array of ploidy for parent 1*/	
}

parameters{
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
}

model{

	for(i in 1:L){
		for(j in 1:N0){
			/* increment likelihood */
			target += calc_lik(GL00[j,i], GL01[j,i], GL02[j,i], P0[i], ploidy0[j,i]);
		}
		for(j in 1:N1){
			/* increment likelihood */
			target += calc_lik(GL10[j,i], GL11[j,i], GL12[j,i], P1[i], ploidy1[j,i]);
		}
	}
	for(i in 1:L){
		/* increment priors on allele frequencies */
		target += beta_lpdf(P0[i] | 0.5, 0.5);
		target += beta_lpdf(P1[i] | 0.5, 0.5);
	}
}

