
data{
	int L; /* # of loci */
	int N; /* # of organisms */
	int J; /* # of populations */
    int<lower=0> popid[N]; /* population ID index */
    real geo[J]; /* 1D array of geographic locations */
	int<lower=0, upper=2> G[N, L]; /* 2D array of G */
}

parameters{
	real<lower=0, upper=1> P[J, L]; /* 2D array of population allele frequencies */
    vector[L] slope; /* vector of slopes */
    vector[L] cent; /* vector of intercepts */
    real mu; /* mean for widths, w */
    real<lower=0> sigma; /* sigma for slopes */
    real<lower=0> serr; /* residual SD */
}

transformed parameters{
    vector[L] w; /* vector of cline widths */
    real sdg; /* default prior on center */
    /* compute cline width slope on the logit scale */
    for(i in 1:L){
        w[i] = 1/(0.25 * slope[i]);
    }
    sdg = 2 * sd(geo);
    
}    

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood for genotype*/
			target += binomial_lpmf(G[j,i] | 2, P[popid[j],i]);
		}
	}
	for(i in 1:L){
	    for(j in 1:J){
    		/* increment linear models */
    		target += normal_lpdf(logit(P[j,i]) | cent[i] + slope[i] * geo[j], serr);
    	}	
	}
	/* increment remaining priors */
	for(i in 1:L){
	    target += normal_lpdf(slope[i] | mu, sigma);
	    target += normal_lpdf(cent[i] | 0, sdg);
	}
	target += normal_lpdf(mu | 0, 5);
	target += gamma_lpdf(sigma | 0.1, 0.01);
	target += gamma_lpdf(serr | 0.1, 0.01);    
}

