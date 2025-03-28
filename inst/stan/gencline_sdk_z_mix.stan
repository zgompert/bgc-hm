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
		} else if(pl==1){ /* haploid locus x ind */
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
	array[N, L] real<lower=0, upper=2> Z; /* 2D array of ancestry*/
	array[N, L] real<lower=0, upper=2> ploidy; /* 2D array of ploidy */	
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
	real<lower=0> sc; /* sigma for center*/
	real<lower=0> sv; /* sigma for v*/
}

parameters{
	vector<lower=0.001,upper=0.999>[L] center; /* cline center parameter */
	vector<lower=0.1,upper=10>[L] v; /* cline width parameter */
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
			target += calc_lik(Z[j,i], H[j], v[i], u[i], ploidy[j,i]);
		}
	}
	for(i in 1:(L)){
		/* increment prior on v and u */
		target += normal_lpdf(logit(center[i]) | 0, sc); /* make sure 0 works for both */
		target += normal_lpdf(log10(v[i]) | 0, sv);

	}
}

