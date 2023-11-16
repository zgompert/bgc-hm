functions {
	/* computes likelihood weighted by genotype likelihoods */
	real calc_lik(real gl0, real gl1, real gl2, real p){
		real prob;
        	prob = log(gl0) + log(1-p) + log(1-p); 
       		prob = log_sum_exp(prob, log(gl1) + log(2) + log(p) + log(1-p));                  
        	prob = log_sum_exp(prob, log(gl2) + log(p) + log(p));
		return prob;
	}
}
data{
	int L; /* # of loci */
	int N0; /* # of parent 0 organisms */
	int N1; /* # of parent 1 organisms */	
	real<lower=0, upper=2> GL00[N0, L]; /* 2D array of genlik 0 for parent 0 */
	real<lower=0, upper=2> GL10[N1, L]; /* 2D array of genlik 0 for parent 1 */
	real<lower=0, upper=2> GL01[N0, L]; /* 2D array of genlik 1 for parent 0 */
	real<lower=0, upper=2> GL11[N1, L]; /* 2D array of genlik 1 for parent 1 */
	real<lower=0, upper=2> GL02[N0, L]; /* 2D array of genlik 2 for parent 0 */
	real<lower=0, upper=2> GL12[N1, L]; /* 2D array of genlik 2 for parent 1 */	
}

parameters{
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
}

model{

	for(i in 1:L){
		for(j in 1:N0){
			/* increment likelihood */
			target += calc_lik(GL00[j,i], GL01[j,i], GL02[j,i], P0[i]);
		}
		for(j in 1:N1){
			/* increment likelihood */
			target += calc_lik(GL10[j,i], GL11[j,i], GL12[j,i], P1[i]);
		}
	}
	for(i in 1:L){
		/* increment priors on allele frequencies */
		target += beta_lpdf(P0[i] | 0.5, 0.5);
		target += beta_lpdf(P1[i] | 0.5, 0.5);
	}
}

