#' The 'analysis2Dmito packge'
#'
#' @description analysis2Dmito fits a Bayesian hierarchical linear mixture model, to
#' OXPHOS protein abundance data. The model is proposed in the paper
#' `Bayesian classification of OXPHOS skeletal myofibres` by Childs et al.. The functions
#' provided here are specific for this purpose. Unforetunately, due to
#' lack of consistency in the structure of related data, tthe user must
#' first wrangle their data into a specified format. A detailed
#' description of the this and how to use the important fucntions within
#' the package can be found on the package
#' [GitHub page](https://github.com/jordanbchilds/analysis2Dmito).
#'
#' @docType package
#' @name analysis2Dmito
#' @aliases analysis2Dmito
#' @useDynLib analysis2Dmito
#' @import methods
#' @import Rcpp
#' @import RcppParallel
#' @import rstan
#' @import rstantools
#' @import data.table
#' @import plyr
#' @import readr
#' @import tidyr
#'
"_PACKAGE"
NULL
