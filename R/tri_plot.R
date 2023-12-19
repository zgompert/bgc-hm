#' Plots interpopulation ancestry (Q10) as a function of hybrid index
#'
#' Creates a triangle plot of hybrid index versus interpopulation ancestry.
#' @param hi a vector of hybrid index estimates (from est_h or est_Q)
#' @param Q10 a vector of interpopulation ancestry estimates (from est_Q)
#' @param pdf a logical specifying whether results should be output to a pdf file; if false the plot is sent to the default graphics device.
#' @param outf a character string specifying the name of the output file if 'pdf=TRUE' [default = tri_plot.pdf].
#' @param ... additional arguments for plotting, see options in par and plot.
#'
#' @details This function generates a scatterplot of interpopulation (a.k.a. interclass or interspecies) ancestry as a function of hybrid index. In other words, this shows the proportion of the genome where each putative hybrid inherited a gene copy from both parents, versus the proportion of the genome inherited from parent 1. Theoretical maxima for interpopulation ancestry given a value of hybrid index are shown as a triangle. Individuals with maximal values of interpopulation ancestry given their hybrid index have one or more non-hybrid parents, meaning they are F1s (hybrid index = 0.5 and interpopulation ancestry = 1) or backcrosses (other cases of maximal interpopulation ancestry). Of course, uncertainty in these admixture parameters can pull point estimates away from these theoretical expectations. 
#'
#' @return A plot is produced, but there is no return value.
#'
#' @export
tri_plot<-function (hi=NULL, Q10=NULL, pdf=TRUE, outf="tri_plot.pdf",...){ 
    if (is.null(hi) == TRUE | is.null(Q10) == TRUE) 
        stop("error, input data are not provided")
    if (pdf == TRUE) 
        pdf(file = paste(out.file))
    plot(hi, Q10, xlab = "Hybrid index", ylab = "Interpopulation ancestry", 
        xlim = c(0, 1), ylim = c(0, 1),...)
    lines(c(0, 0.5), c(0, 1))
    lines(c(0.5, 1), c(1, 0))
    if (pdf == TRUE) 
        dev.off()
}

