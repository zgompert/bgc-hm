functions {
	real calc_lik(real g, real p0, real p1, real q11, real q10, real q00){
		real prob;
		if(g==0)
			prob = log(q11 * (1-p1) * (1-p1) + q00 * (1-p0) * (1-p0) + q10 * (1-p1) * (1-p0));
		else if (g==1)
			
			prob = log(q11 * 2 * (1-p0) * p0 + q00 * 2 * (1-p1) * p1 + q10 *  (1-p0) * p1 + q10 * p0 * (1-p1));
		else
			prob = log(q11 * p1 * p1 + q00 * p0 * p0 + q10 * p1 * p0);
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	int<lower=0, upper=2> G[N, L]; /* 2D array of G */
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
}

parameters{
	simplex[3] Q[N]; /* array of simplexes for ancestry components */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(G[j,i], P0[i], P1[i], Q[j][1], Q[j][2], Q[j][3]);
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

