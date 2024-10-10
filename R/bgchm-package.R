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
#' Gompert Z, DeRaad D, Buerkle CA. A next generation of hierarchical Bayesian analyses of hybrid zones enables model-based quantification of variation in introgression in R. bioRxiv 2024.03.29.587395.
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.7. https://mc-stan.org
#'
NULL
