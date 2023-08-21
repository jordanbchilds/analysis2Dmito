#' @title Get Example Dataset
#'
#' @description
#' Returns an example data set. The dataset consists of protein expression levels
#' collected via IMC, for VDAC, NDUFB8, CYB and MTCO1, on 12 patients.
#'
#' @return A [data.frame] object.
#'
#' @examples
#' exampleData = get_exampleData()
#'
#' @description
#' The data is returned as a [data.frame] object in the correct form required for
#' use in the [anaylsis2Dmito::getData_mats] function.
#'
#' @details
#' Pre-processing of the data here is only to remove values which are less than
#' or equal to zero. This so that the data can be logged transformed without
#' problem. Assuming the original dataset will not be updated this results in
#' the removal of two fibres.
#'
#' @author Amy Vincent et al.
#'
#' @importFrom readr read_delim
#' @importFrom tidyr pivot_longer
#' @importFrom plyr match_df
#'
#' @export
get_exampleData = function() {
  urlfile = "https://raw.githubusercontent.com/jordanbchilds/AV_mitocyto/Data_prepped,.csv"
  rawData = readr::read_delim(url(urlfile))

  rawData = rawData[,c("ID", "patient_id", mitochan, channels)]
  colnames(rawData) = c("fibreID", "sampleID", mitochan, channels)

  data = tidyr::pivot_longer(rawData, cols=c("VDAC", "NDUFB8", "CYB", "MTCO1"), names_to="channel")

  data_df= as.data.frame(data)

  return( data_df )
}











