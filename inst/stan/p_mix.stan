
data{
	int L; /* # of loci */
	int N0; /* # of parent 0 organisms */
	int N1; /* # of parent 1 organisms */	
	int<lower=0, upper=2> G0[N0, L]; /* 2D array of G for parent 0 */
	int<lower=0, upper=2> G1[N1, L]; /* 2D array of G for parent 1 */
	int<lower=0, upper=2> ploidy0[N0, L]; /* 2D array of ploidy for parent 0 */
	int<lower=0, upper=2> ploidy1[N1, L]; /* 2D array of ploidy for parent 1*/	
}

parameters{
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
}

model{

	for(i in 1:L){
		for(j in 1:N0){
			/* increment likelihood */
			target += binomial_lpmf(G0[j,i] | ploidy0[j,i], P0[i]);
		}
		for(j in 1:N1){
			/* increment likelihood */
			target += binomial_lpmf(G1[j,i] | ploidy1[j,i], P1[i]);
		}
	}
	for(i in 1:L){
		/* increment priors on allele frequencies */
		target += beta_lpdf(P0[i] | 0.5, 0.5);
		target += beta_lpdf(P1[i] | 0.5, 0.5);
	}
}

