
data{
	int J; /* # of populations */
	real geo[J]; /* 1D array of geographic locations */
	real Y[J]; /* 1D array of logit population allele frequencies */	
	real lb; /* lower bound of logit p to include */
	real ub; /* lower bound of logit p to include */
}

transformed data{
    real sdg; /* default prior on center */
    sdg = 3 * sd(geo);
}

parameters{
	
    real slope; /* slope */
    real cent; /* intercept */
    real<lower=0> serr; /* residual SD */
}

transformed parameters{
    real w; /* cline width */
    /* compute cline width slope on the logit scale */
    w = 1/(0.25 * slope);
    
}    

model{
	for(j in 1:J){
    	/* increment likelihood, only if logit p is in range */
		if((Y[j] > lb) && (Y[j] < ub))
			target += normal_lpdf(Y[j] | cent + slope * geo[j], serr);
    		
	}
	/* increment priors */
	target += normal_lpdf(slope | 0, 10);
	target += normal_lpdf(cent | 0, sdg);
	target += gamma_lpdf(serr | 0.1, 0.01);    
}

