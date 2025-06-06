functions {
	real calc_lik(real g, real p0, real p1, real h){
		real prob;
		if(g==0)
			prob = log(h * (1-p1) + (1-h) * (1-p0)) + log(h * (1-p1) + (1-h) * (1-p0));
		else if (g==1)
			prob = log(2) + log(h * (1-p1) + (1-h) * (1-p0)) + log(h * p1 + (1-h) * p0);
		else
			prob = log(h * p1 + (1-h) * p0) + log(h * p1 + (1-h) * p0);
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	array[N, L] int<lower=0, upper=2> G; /* 2D array of G */
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
			target += calc_lik(G[j,i], P0[i], P1[i], H[j]);
		}
	}
	for(j in 1:N){
		/* increment prior on hybrid index */
		target += beta_lpdf(H[j] | 0.5, 0.5);

	}
}

