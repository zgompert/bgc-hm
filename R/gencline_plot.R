#' Plots genomic clines for a set of loci
#'
#' Plots a set of genomic clines.
#' @param center vector of cline centers (from est_gencline.R)
#' @param v vector of cline gradients (from est_gencline.R) 
#' @param pdf a logical specifying whether results should be output to a pdf file
#' @param outf a character string specifying the name of the output file if 'pdf=TRUE'
#' @param ... additional arguments for plotting
#'
#' @return A plot is produced, but there is no return value
#'
#' @export
gencline_plot<-function (center=NULL, v=NULL, pdf=TRUE, outf="tri_plot.pdf", ...){ 
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
    for(i in 1:length(v)){
	    phi<-(h^v[i])/(h^v[i] + (1-h)^v[i] * exp(u[i]))
	    lines(h,phi, ...)
    }
    abline(a=0,b=1,lty=2)
    if (pdf == TRUE) 
        dev.off()
}

