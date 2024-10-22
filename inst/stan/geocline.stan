
data{
	int L; /* # of loci */
	int J; /* # of populations */
	array[J] real geo; /* 1D array of geographic locations */
	array[J, L] real Y; /* 2D array of logit population allele frequencies */	
	real lb; /* lower bound of logit p to include */
	real ub; /* lower bound of logit p to include */
	real ga; /* alpha parameter for gamma priors */
	real gb; /* beta parameter for gamma priors */
}

transformed data{
    real sdg; /* default prior on center */
    real sds; /* deafult prior on slope */
    sdg = 3 * sd(geo); /* based on empirical SD of geo */
    sds = 1.5 * 4/sd(geo);
}

parameters{
	
    vector[L] slope; /* vector of slopes */
    vector[L] cent; /* vector of intercepts */
    real mu; /* mean for widths, w */
    real<lower=0> sigma; /* sigma for slopes */
    real<lower=0> serr; /* residual SD */
}

transformed parameters{
    vector[L] w; /* vector of cline widths */
    /* compute cline width slope on the logit scale */
    for(i in 1:L){
        w[i] = 1/(0.25 * slope[i]);
    }
    
}    

model{
	for(i in 1:L){
	    for(j in 1:J){
    		/* increment likelihood, only if logit p is in range */
		if((Y[j,i] > lb) && (Y[j,i] < ub))
	    		target += normal_lpdf(Y[j,i] | cent[i] + slope[i] * geo[j], serr);
    	    }	
	}
	/* increment priors */
	for(i in 1:L){
	    target += normal_lpdf(slope[i] | mu, sigma);
	    target += normal_lpdf(cent[i] | 0, sdg);
	}
	target += normal_lpdf(mu | 0, sds);
	target += gamma_lpdf(sigma | ga, gb);
	target += gamma_lpdf(serr | ga, gb);    
}

