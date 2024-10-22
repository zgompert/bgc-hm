functions {
	real calc_lik(real gl0, real gl1, real gl2, real p0, real p1, real q11, real q10, real q00){
		real prob;
		prob = log(gl0) + log(q11 * (1-p1) * (1-p1) + q00 * (1-p0) * (1-p0) + q10 * (1-p1) * (1-p0));
		prob = log_sum_exp(prob, log(gl1) + log(q11 * 2 * (1-p0) * p0 + q00 * 2 * (1-p1) * p1 + q10 *  (1-p0) * p1 + q10 * p0 * (1-p1)));
	    prob = log_sum_exp(prob, log(gl2) + log(q11 * p1 * p1 + q00 * p0 * p0 + q10 * p1 * p0));
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	array[N, L] real<lower=0, upper=2> GL0; /* 2D array of genlik 0*/
	array[N, L] real<lower=0, upper=2> GL1; /* 2D array of genlik 1*/
	array[N, L] real<lower=0, upper=2> GL2; /* 2D array of genlik 2*/
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
}

parameters{
	array[N] simplex[3] Q; /* array of simplexes for ancestry components */
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(GL0[j,i], GL1[j,i], GL2[j,i], P0[i], P1[i], Q[j][1], Q[j][2], Q[j][3]);
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

