#' Plots genomic clines for a set of loci
#'
#' Plots a set of genomic clines.
#' @param center vector of cline centers (from est_gencline)
#' @param v vector of cline gradients (from est_gencline) 
#' @param pdf a logical specifying whether results should be output to a pdf file; if false the plot is sent to the default graphics device. 
#' @param outf a character string specifying the name of the output file if 'pdf=TRUE' [default = cline_plot.pdf].
#' @param cvec vector of colors to use for clines (one entry per locus)
#' @param xvec vector of line widths (relative to 1) to use for clines (one entry per locus)
#' @param ... additional arguments for plotting, see options in par and plot.
#'
#' @details This function plots genomic clines for set of loci, that is the probability of local ancestry from parental population 1 at a locus given hybrid index (the overall proportion of an individual's genome inherited from population 1). The clines for all loci are shown on a single plot, with one line per locus. A 1:1 dashed line denotes the null expected ancestry probability if all loci exhibit dynamics precisely equal to the genome-wide average intogression.
#'
#' @return A plot is produced, but there is no return value.
#'
#' @export
gencline_plot<-function (center=NULL, v=NULL, pdf=TRUE, outf="cline_plot.pdf", cvec=NULL, xvec=NULL,
			 ...){ 
    if (is.null(center) == TRUE | is.null(v) == TRUE) 
        stop("error, input data are not provided")
    if (pdf == TRUE) 
        pdf(file = paste(out.file))

    ## create vector for hybrid indexes
    h<-seq(0,1,0.01)
    ## convert center to u
    u<-log(center/(1-center)) * v
	    
    plot(h, h, type='n', xlab = "Hybrid index", ylab = "Ancestry probability", 
        xlim = c(0, 1), ylim = c(0, 1), ...)
    if(is.null(cvec) & is.null(xvec)){
	    for(i in 1:length(v)){
		 phi<-(h^v[i])/(h^v[i] + (1-h)^v[i] * exp(u[i]))
	   	 lines(h,phi, ...)
    	    }
    } else if (is.null(cvec) & is.null(xvec)==FALSE){
	    cvec<-rep("darkgray",length(v))
	    for(i in 1:length(v)){
		 phi<-(h^v[i])/(h^v[i] + (1-h)^v[i] * exp(u[i]))
	   	 lines(h,phi,col=cvec[i],lwd=xvec[i])
    	    }
    } else if (is.null(cvec)==FALSE & is.null(xvec)){
	    xvec<-rep(1,length(v))
	    for(i in 1:length(v)){
		 phi<-(h^v[i])/(h^v[i] + (1-h)^v[i] * exp(u[i]))
	   	 lines(h,phi,col=cvec[i],lwd=xvec[i])
    	    }
    } else{
            for(i in 1:length(v)){
                 phi<-(h^v[i])/(h^v[i] + (1-h)^v[i] * exp(u[i]))
                 lines(h,phi,col=cvec[i],lwd=xvec[i])
            }   
    }

    abline(a=0,b=1,lty=2)
    if (pdf == TRUE) 
        dev.off()
}

