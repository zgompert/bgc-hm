functions {
	real calc_lik(real gl0, real gl1, real gl2, real p0, real p1, real h, real pl){
		real prob;
		if(pl==2){/* diploid locus x ind */
    		prob = log(gl0) + log(h * (1-p1) + (1-h) * (1-p0)) + log(h * (1-p1) + (1-h) * (1-p0));
	    	prob = log_sum_exp(prob, log(gl1) + log(2) + log(h * (1-p1) + (1-h) * (1-p0)) + log(h * p1 + (1-h) * p0));
	    	prob = log_sum_exp(prob, log(h * p1 + (1-h) * p0) + log(h * p1 + (1-h) * p0));
	    } else if(pl==1){ /* haploid locus x ind */
	        prob = log(gl0) + log(h * (1-p1) + (1-h) * (1-p0));
	    	prob = log_sum_exp(prob, log(gl1) + log(h * p1 + (1-h) * p0));
	    } else{
	    	prob=0;
	    }
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	array[N, L] real<lower=0, upper=2> GL0; /* 2D array of genlik 0*/
	array[N, L] real<lower=0, upper=2> GL1; /* 2D array of genlik 1*/
	array[N, L] real<lower=0, upper=2> GL2; /* 2D array of genlik 2*/
	array[N, L] real<lower=0, upper=2> ploidy; /* 2D array of ploidy */
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
}

parameters{
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(GL0[j,i], GL1[j,i], GL2[j,i], P0[i], P1[i], H[j], ploidy[j,i]);
		}
	}
	for(j in 1:N){
		/* increment prior on hybrid index */
		target += beta_lpdf(H[j] | 0.5, 0.5);

	}
}

