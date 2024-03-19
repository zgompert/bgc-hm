#' Simulated hybrid zone data with known genotypes 
#'
#' @description This data set was created from a hybrid zone simulation with dfuse. This data set comprises 51 ancestry-informative loci, 100 putative hybrids and 50 individuals representative of each parental species. All loci are diploid. The data were simulated with dfuse assuming 110 demes, m = 0.1 between neighboring demes, and 10 impacting hybrid fitness via underdominance. The data set comprises 3 matrixes, GenHybrids, GenP0, and GenP1, which contain the genotypic data (biallelic SNPs) for the hybrids, parental individuals respectively.
#'
#' @format GenHybrids, GenP0 and GenP1 each constitute a \code{matrix} with 51 columns (one per locus) and N rows (one per individual). N = 100 for the hybrids and 50 for each parental population. Each element in each matrix contains the genotype for an individual, coded as 0, 1 or 2, where 0 and 2 denote alternative homozygotes and 1 denotes heterozygotes.
#'
#' @source This is simulated data that can be used for example or test analyses.
#'
#' @export 
