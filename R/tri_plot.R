#' Plots interpopulation ancestry (Q10) as a function of hybrid index
#'
#' Uses Hamiltonian Monte Carlo (HMC) for Bayesian inference of ancestry classes from genetic data. Ancestry classes denote the proportion of an individual's genome where both gene copies come from source 1 (Q11), both gene copies come from source 0 (Q00), or where one gene copy comes from source 1 and one from source 0.
#' @param hi a vector of hybrid index estimates (from est_h or est_Q)
#' @param Q10 a vector of interpopulation ancestry estimates (from est_Q)
#' @param pdf a logical specifying whether results should be output to a pdf file
#' @param outf a character string specifying the name of the output file if 'pdf=TRUE'
#'
#' @return A plot is produced, but there is no return value
#'
#' @export
tri_plot<-function (hi=NULL, Q10=NULL, pdf=TRUE, outf="tri_plot.pdf"){ 
    if (is.null(hi) == TRUE | is.null(Q10) == TRUE) 
        stop("error, input data are not provided")
    if (pdf == TRUE) 
        pdf(file = paste(out.file))
    plot(hi, Q10, xlab = "Hybrid index", ylab = "Interpopulation ancestry", 
        xlim = c(0, 1), ylim = c(0, 1))
    lines(c(0, 0.5), c(0, 1))
    lines(c(0.5, 1), c(1, 0))
    if (pdf == TRUE) 
        dev.off()
}
