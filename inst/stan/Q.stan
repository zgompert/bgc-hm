functions {
	real calc_lik(real g, real p0, real p1, real q11, q12, q22){
		real prob;
		if(g==0)
			prob = log(q11 * (1-p1) * (1-p1) + q22 * (1-p0) * (1-p0) + q12 * 2 * (1-p1) * (1-p0));
		else if (g==1)
			
			prog = log(q11 * 2 * (1-p0) * p0 + q22 * 2 * (1-p1) * p1 + q12 *  (1-p0) * p1 + q12 * p0 * (1-p1));
		else
			prob = log(q11 * p1 * p1 + q22 * p0 * p0 + q12 * 2 * p1 * p0);
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	real<lower=0, upper=2> G[N, L]; /* matrix of G */
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 0 allele frequencies */
}

parameters{
	matrix<lower=0, upper=1>[N, 3] Q; /* matrix of ancestry components */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(G[j,i], P0[i], P1[i], Q[j,1], Q[j,2], Q[j,3]);
		}
	}
	for(j in 1:N){
		/* increment prior on Q */
		target += beta_lpdf(Q[j,] | rep_vector(0.5, 3));

	}
}

