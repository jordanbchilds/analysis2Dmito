#' @title Example QIF dataset
#'
#' @description
#' This dataset has been taken from \href{https://www.medrxiv.org/content/10.1101/2022.04.04.22272484v1.full}{"Estimating the allele frequency threshold of
#' the pathogenic mtDNA variant m.3243A>G tolerated by human myofibres"} by Ahmed
#' et al. (2022).
#'
#' @details
#' The data is read from the \href{https://github.com/CnrLwlss/Ahmed_2022}{GitHub repository}
#' associated with the paper and so any changes made to the data there would
#' remain here. Preprocessing of the data here is only to remove values which
#' are less than or equal to zero. This so that the data can be logged
#' transformed without problem. Assuming the original dataset will not be
#' updated this results in the removal of two fibres.
#'
#' @author Syeda T Ahmed, Robert W Taylor, Doug M Turnbull, Conor Lawless, Sarah T Pickett
#'
#' @importFrom readr read_delim
#' @importFrom tidyr pivot_longer
#' @importFrom plyr match_df
"exampleData"
