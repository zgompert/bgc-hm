functions {
	real calc_phi(real h, real vv, real uu){
		real phi;
		phi = (h^vv)/((h^vv)+((1-h)^vv)*exp(uu));
		return phi;
	}
	real calc_lik(real z, real h, real vv, real uu, real pl){
		real prob;
		real phi;
		phi = calc_phi(h, vv, uu);
		if(pl==2){/* diploid locus x ind */
        		if(z==0)
            			prob = log(1-phi) + log(1-phi);
        		else if (z==1)
        	    		prob = log(2) + log(phi) + log(1-phi);
       		 	else
            			prob = log(phi) + log(phi);
		} else if (pl==1){ /* haploid locus x ind */
			if(z==0)
            			prob = log(1-phi);
       		 	else
            			prob = log(phi);
            	} else{
            		prob=0;
            	}		
		return prob;
	}

}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	real<lower=0, upper=2> Z[N]; /* 1D array of ancestry*/
	real<lower=0, upper=2> ploidy[N]; /* 1D array of ploidy */
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
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
		target += calc_lik(Z[j], H[j], v, u, ploidy[j]);
	}

	/* increment prior on v and u */
	target += normal_lpdf(logit(center) | 0, sc); 
	target += normal_lpdf(log10(v) | 0, sv);
}

