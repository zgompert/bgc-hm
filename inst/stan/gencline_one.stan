functions {
	real calc_phi(real h, real vv, real uu){
		real phi;
		phi = (h^vv)/((h^vv)+((1-h)^vv)*exp(uu));
		return phi;
	}
	real calc_lik(real g, real p0, real p1, real h, real vv, real uu){
		real prob;
		real phi;
		phi = calc_phi(h, vv, uu);
        if(g==0)
            prob = log(phi * (1-p1) + (1-phi) * (1-p0)) + log(phi * (1-p1) + (1-phi) * (1-p0));
        else if (g==1) 
            prob = log(phi * (1-p1) + (1-phi) * (1-p0)) + log(phi * p1 + (1-phi) * p0);
        else    
            prob = log(phi * p1 + (1-phi) * p0) + log(phi * p1 + (1-phi) * p0);

		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	array[N] real<lower=0, upper=2> G; /* 1D array of G*/
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
	real<lower=0, upper=1> P0; /* parent 0 allele frequencies */
	real<lower=0, upper=1> P1; /* parent 1 allele frequencies */
	real<lower=0> sc; /* sigma for center*/
	real<lower=0> sv; /* sigma for v*/
}

parameters{
	real<lower=0.001,upper=0.999> center; /* cline center parameter */
	real<lower=0.1,upper=10> v; /* cline width parameter */
}

transformed parameters{
	real<lower=-100, upper=100> u; /* cline u parameter */
	u = logit(center) * v;
}

model{

	for(j in 1:N){
		/* increment likelihood */
		target += calc_lik(G[j], P0, P1, H[j], v, u);
	}

	/* increment prior on v and u */
	target += normal_lpdf(logit(center) | 0, sc); 
	target += normal_lpdf(log10(v) | 0, sv);
}

