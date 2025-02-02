#' The 'decemedip' package.
#'
#' @description The R package decemedip is a novel computational paradigm developed
#' for inferring the relative abundances of cell types and tissues from tissue bulk
#' or circulating cell-free DNA (cfDNA) measure by methylated DNA immunoprecipitation
#' sequencing (MeDIP-Seq). This paradigm allows using reference data from other
#' technologies such as microarray or WGBS.
#'
#' @name decemedip-package
#' @aliases decemedip-package
#' @useDynLib decemedip, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.6. https://mc-stan.org
#'
"_PACKAGE"
