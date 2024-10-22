functions {
	real calc_phi(real h, real vv, real uu){
		real phi;
		phi = (h^vv)/((h^vv)+((1-h)^vv)*exp(uu));
		return phi;
	}
	real calc_lik(real z, real h, real vv, real uu){
		real prob;
		real phi;
		phi = calc_phi(h, vv, uu);
        if(z==0)
            prob = log(1-phi) + log(1-phi);
        else if (z==1)
            prob = log(2) + log(phi) + log(1-phi);
        else
            prob = log(phi) + log(phi);
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	array[N, L] real<lower=0, upper=2> Z; /* 2D array of ancestry*/
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
    real<lower=0> sd0; /* SD for SD of normal prior on cline distributions*/
    real<lower=0> mu0; /* SD for mean of normal prior on cline distributions*/
}

parameters{
	vector<lower=0.001,upper=0.999>[L] center; /* cline center parameter */
	vector<lower=0.1,upper=10>[L] v; /* cline width parameter */
	real<lower=0> sc; /* sigma for center*/
	real<lower=0> sv; /* sigma for v*/
    real<lower=0> muc; /* mean for center*/
    real<lower=0> muv; /* mean for v*/

}

transformed parameters{
	vector<lower=-100, upper=100>[L] u; /* cline u parameter */
	for(i in 1:L){
		u[i] = logit(center[i]) * v[i];
	}
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(Z[j,i], H[j], v[i], u[i]);
		}
	}
	for(i in 1:(L)){
		/* increment prior on v and u */
		target += normal_lpdf(logit(center[i]) | muc, sc); 
		target += normal_lpdf(log10(v[i]) | muv, sv);

	}
	/* increment prior on sc and sv */
	target += normal_lpdf(sc | 0, sd0);
	target += normal_lpdf(sv | 0, sd0);
        /* increment prior on muc and muv */
        target += normal_lpdf(muc | 0, mu0);
        target += normal_lpdf(muv | 0, mu0);

}

generated quantities{
	real<lower=0> sdc; /* actual sd of cline params */
	real<lower=0> sdv; /* actual sd of cline params */

	sdc = sd(logit(center));
	sdv = sd(log10(v));
}

