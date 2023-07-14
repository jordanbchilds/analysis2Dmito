#' @title Get Example Dataset
#'
#' @description
#' Returns an example data set. The dataset has been taken from
#' \href{https://www.medrxiv.org/content/10.1101/2022.04.04.22272484v1.full}{"Estimating t
#' he allele frequency threshold of the pathogenic mtDNA variant m.3243A>G
#' tolerated by human myofibres"} by Ahmed et al. (2022). It is read directly
#' the Github repository associated with the paper.
#'
#' @return A [data.frame] object.
#'
#' @examples
#' exampleData = get_exampleData()
#'
#' @description
#' This dataset has been taken from \href{https://www.medrxiv.org/content/10.1101/2022.04.04.22272484v1.full}{"Estimating the allele frequency threshold of
#' the pathogenic mtDNA variant m.3243A>G tolerated by human myofibres"} by Ahmed
#' et al. (2022).
#'
#' @details
#' Preprocessing of the data here is only to remove values which are less than
#' or equal to zero. This so that the data can be logged transformed without
#' problem. Assuming the original dataset will not be updated this results in
#' the removal of two fibres.
#'
#' @author Syeda T Ahmed, Robert W Taylor, Doug M Turnbull, Conor Lawless, Sarah T Pickett
#'
#' @importFrom readr read_delim
#' @importFrom tidyr pivot_longer
#' @importFrom plyr match_df
#'
#' @export
get_exampleData = function() {
  urlfile = "https://raw.githubusercontent.com/CnrLwlss/Ahmed_2022/master/rawdat_filtered_a.csv"
  rawData = readr::read_delim(url(urlfile))

  rawData_sub = rawData[, c("caseno", "Fibre", "raw_porin", "raw_CI", "raw_CIV", "controls")]
  longForm_df = tidyr::pivot_longer(rawData_sub,
                                    cols=c("raw_porin", "raw_CI", "raw_CIV"),
                                    names_to="channels")
  colnames(longForm_df) = c("sampleID", "fibreID", "sbj_type", "channel", "value")

  longForm_df = as.data.frame(longForm_df)

  longForm_df$unique_id = paste(longForm_df$sampleID, longForm_df$fibreID, sep="_")
  zeroFibs = longForm_df[longForm_df$value == 0.0, "unique_id"]
  prepped_data = longForm_df[!(longForm_df$unique_id %in% zeroFibs), ]
  prepped_data$unique_id = NULL

  return(prepped_data)
}











