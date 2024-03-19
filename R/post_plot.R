#' Plots bivariate posterior distribution of hierarchical cline parameters 
#'
#' Creates bivariate plots of combinations of cline standard deviations and means for one or multiple sets of loci. 
#' @param objs the output from est_gencline, or a list of such objects 
#' @param param1 a string specifying the genomic cline parameter (muc, muv, sc or sv) for the x-axis on the plot. 
#' @param param2 a string specifying the genomic cline parameter (muc, muv, sc or sv) for the y-axis on the plot. 
#' @param probs vector specifying which confidence ellipses to plot (e.g., 0.95 for an ellipse that encloses 95% of the posterior).
#' @param colors single color or vector of colors for plotting, if more than one color is supplied different colors will be used for different posteriors (in the order specified by objs).
#' @param addPoints Boolean indicating whether points denoting samples from the posterior should be added to the plot (default = TRUE).
#' @param palpha transparency (alpha) value for the color of points, must be between 0 and 1, smaller values make points more transparent.
#' @param pdf a logical specifying whether results should be output to a pdf file; if false the plot is sent to the default graphics device.
#' @param outf a character string specifying the name of the output file if 'pdf=TRUE' [default = tri_plot.pdf].
#' @param ... additional arguments for plotting, see options in par and plot.
#'
#' @details This function produces a plot illustrating the bivariate posterior probability distribution for the cline mean or standard deviation parameters for one or more sets of genetic loci. Bivariate plots can be produced for any combination of higher-level cline parameters, that is the mean center (muc), mean slope (muv), standard deviation of centers (sc) or standard deviation of slopes (sv) (these are the names of the parameters in the HMC objects and should be used for param1 and param2). Plots include confidence ellipses cpaturing user specified proportions of the posterior (set with probs) and optionally points representing the samples from the posterior. Ellipses are produced by approximating the bivariate posterior with a bivariate normal distribution (these are meant to serve an illustrative purpose only). Posteriors from multiple sets of loci (multiple outputs from the est_gencline function) can be placed on the same plot. To plot a single posterior, pass the entire object (a list) produced by est_gencline to the function for the objs argument. To plot multiple, pass a list of these objects (a list of lists) to the objs function.
#'
#' @details Users can specify colors for plotting posteriors. When more than one color is provided, different colors will be used for each object (set of loci). The transparency of the points (if included) can be adjusted with the palpha parameter. Transperancy of the ellipses is set automatically by setting the alpha (transparency) value for each ellipses to 1-probs.
#'
#' @return A plot is produced, but there is no return value.
#'
#' @export
pp_plot<-function(objs=NULL,param1="muc",param2="sdc",probs=c(0.5,0.75,0.95),colors="black",addPoints=TRUE,palpha=0.3,pdf=TRUE, outf="pp_plot.pdf",...){ 
	## extract relevant parameters from the hmc objects
	if(is.list(objs)==FALSE){
		stop("a list object must be supplied to objs")
	} 
	if("gencline_hmc" %in% names(objs)){ ## single list object
		N<-1
		ojbs<-list(objs,NULL) ## just to make the code work in a similar way to the case where there are multiple object
	} else{
		N<-length(objs)
	}
	## set up colors for plotting
	colors<-rep(colors,length.out=N)
	probs<-sort(probs)
	alphas<-seq(.5,1,length.out=length(probs))
	clist<-vector("list",N)
	for(k in 1:N){
		clist[[k]]<-rep(NA,length(probs))
		for(j in length(probs)){
			rgb <- farver::decode_colour(colors[k], alpha = TRUE)
    			rgb[!is.na(alpha), 4]<-1-probs[j]
    			clist[[k]][j]<-farver::encode_colour(rgb, rgb[, 4])
		}
		rgb <- farver::decode_colour(colors[k], alpha = TRUE)
    		rgb[!is.na(alpha), 4]<-palpha
    		colors[[k]]<-farver::encode_colour(rgb, rgb[, 4])
	}

	## create vectors to store HMC samples
	p1<-vector("list",N)
	p2<-vector("list",N)
	for(k in 1:N){
		p1[[k]]<-rstan::extract(objs[[k]]$gencline_hmc,param1)[[1]]
		p2[[k]]<-rstan::extract(objs[[k]]$gencline_hmc,param2)[[1]]
	}
	## compute plot bounds
	xbnds<-c(min(unlist(p1)),max(unlist(p1)))
	ybnds<-c(min(unlist(p2)),max(unlist(p2)))
	
	if (pdf == TRUE) 
		pdf(file = paste(out.file))
  	plot(0,0,xlim=xbnds,ylim=ybnds,type='n',...)
	for(k in 1:N){
		if(addPoints==TRUE)
			points(p1[[k]],p2[[k]],col=colors[k],...)
		# posterior means and covariance matrixes
		mu <- c(mean(p1[[k]]), mean(p2[[k]]))
		cov_matrix <- matrix(c(var(p1[[k]]), cov(p1[[k]],p2[[k]]), cov(p1[[k]],p2[[k]]), var(p2[[k]])), nrow = 2)
		# compute eigenvalues and vectors, eigenvectors denote axes
		ed <- eigen(cov_matrix)
		eigenvalues <- ed$values
		eigenvectors <- ed$vectors
		for(j in 1:length(probs)){
			# compute scaling factors
			chi2value <- qchisq(probs[j], df = 2)  # df 2 for a bivariate normal distribution
			sf <- sqrt(chi2value * eigenvalues)
			# generate 100 points to create the ellipse
			theta <- seq(0, 2 * pi, length.out = 100)
			eps <- matrix(0, ncol = 2, nrow = length(theta))
			for (i in 1:100) {
				eps[i, ] <- mu + sf[1] * cos(theta[i]) * eigenvectors[, 1] + sf[2] * sin(theta[i]) * eigenvectors[, 2]
			}
			lines(eps,clist[[k]][j],...)
		}
	}
	if (pdf == TRUE) 
		dev.off()
}

