functions {
	real calc_lik(real z, real h, real pl){
		real prob;
		if(pl==2){ /* diplid locus x ind */
        	if(z==0)
            		prob = log(1-h) + log(1-h);
        		else if (z==1)
            		prob = log(2) + log(h) + log(1-h);
        		else
            		prob = log(h) + log(h);
        	} else if(pl==1){ /* haplid locus x ind */
        	        if(z==0)
            		prob = log(1-h);
        		else
            		prob = log(h);
        	} else{
        		prob=0;
        	}
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	int<lower=0, upper=2> Z[N, L]; /* 2D array of ancestry */
	real<lower=0, upper=2> ploidy[N, L]; /* 2D array of ploidy */
}

parameters{
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(Z[j,i], H[j], ploidy[j,i]);
		}
	}
	for(j in 1:N){
		/* increment prior on hybrid index */
		target += beta_lpdf(H[j] | 0.5, 0.5);

	}
}

