functions {
	real calc_phi(real h, real vv, real uu){
		real phi;
		phi = (h^vv)/((h^vv)+((1-h)^vv)*exp(uu));
		return phi;
	}
	/* computes likelihood weighted by genotype likelihoods */
	real calc_lik(real gl0, real gl1, real gl2, real p0, real p1, real h, real vv, real uu, real pl){
		real prob;
		real phi;
		phi = calc_phi(h, vv, uu);
		if(pl==2){/* diploid locus x ind */
            prob = log(gl0) + log(phi * (1-p1) + (1-phi) * (1-p0)) + log(phi * (1-p1) + (1-phi) * (1-p0)); 
            prob += log(gl1) + log(phi * (1-p1) + (1-phi) * (1-p0)) + log(phi * p1 + (1-phi) * p0);                  
            prob += log(gl2) + log(phi * p1 + (1-phi) * p0) + log(phi * p1 + (1-phi) * p0);
		} else { /* haploid locus x ind */
		    prob = log(gl0) + log(phi * (1-p1) + (1-phi) * (1-p0));  
		    /* store one copy alt allele as gl2 */                
            prob += log(gl2) + log(phi * p1 + (1-phi) * p0);
		}
		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	real<lower=0, upper=2> GL0[N, L]; /* 2D array of genlik 0*/
	real<lower=0, upper=2> GL1[N, L]; /* 2D array of genlik 1*/
	real<lower=0, upper=2> GL2[N, L]; /* 2D array of genlik 2*/
	real<lower=0, upper=2> ploidy[N, L]; /* 2D array of ploidy */
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 1 allele frequencies */
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
			target += calc_lik(GL0[j,i], GL1[j,i], GL2[j,i], P0[i], P1[i], H[j], v[i], u[i], ploidy[j,i]);
		}
	}
	for(i in 1:(L)){
		/* increment prior on v and u */
		target += normal_lpdf(logit(center[i]) | 0, sc); /* make sure 0 works for both */
		target += normal_lpdf(log10(v[i]) | 0, sv);

	}
}

