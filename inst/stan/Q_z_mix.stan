functions {
	real calc_lik(real z, real q11, real q10, real q00, real pl){
		real prob;
		if(pl==2){ /* diploid locus x ind*/
		    if(z==0)
		    	prob = log(q00);
		    else if (z==1)
		    	prob = log(q10);
		    else
			    prob = log(q11);
		    return prob;
        } else if(pl==0){ /* haploid locus x ind */
            if(z==0)
                prob = log(q00) + 0.5 * log(q10);
            else
                prob = log(q11) + 0.5 * log(q10);   
        } else{
        	prob=0;
        }
        return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	array[N, L] int<lower=0, upper=2> Z; /* 2D array of ancestry */
    array[N, L] real<lower=0, upper=2> ploidy; /* 2D array of ploidy */    
}

parameters{
	array[N] simplex[3] Q; /* array of simplexes for ancestry components */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(Z[j,i], Q[j][1], Q[j][2], Q[j][3], ploidy[j,i]);
		}
	}
	for(j in 1:N){
		/* increment prior on Q */
		target += dirichlet_lpdf(Q[j] | rep_vector(0.5, 3));

	}
}

generated quantities {
	vector<lower=0, upper=1>[N] H; /* hybrid indexes, computed from Q */

	for(j in 1:N){
		H[j] = Q[j][1] + 0.5 * Q[j][2];
	}
}

