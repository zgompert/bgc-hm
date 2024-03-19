#' The 'bgchm' package.
#'
#' @description This is an R package for Bayesian analyses of population genomic data from hybrid zones, including Bayesian genomic cline analysis, estimation of hybrid indexes and ancestry class proportions, some geographic cline analyses, and accessory plotting functions. This package using Hamiltonian Monte Carlo (HMC) for sampling posterior distributions, with HMC sampling implemented via Stan. 
#'
#' @docType package
#' @name bgchm-package
#' @aliases bgchm
#' @useDynLib bgchm, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Gompert Z, et al. 2024. Bayesian analyses of hybrid zones in R with Hamiltonian Monte Carlo. Manuscript in preparation.
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.7. https://mc-stan.org
#'
NULL
