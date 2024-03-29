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
#' use in the [analysis2Dmito::getData_mats] function.
#'
#'
#' @author Amy Vincent et al.
#'
#' @importFrom readr read_delim
#' @importFrom tidyr pivot_longer
#' @importFrom plyr match_df
#'
#' @export
get_exampleData = function() {
  urlfile = "https://raw.githubusercontent.com/jordanbchilds/AV_mitocyto/main/Data_prepped.csv"
  rawData = readr::read_delim(url(urlfile))

  mitochan = "VDAC"
  channels = c("NDUFB8", "CYB", "MTCO1")

  rawData = rawData[,c("ID", "patient_id", mitochan, channels)]
  colnames(rawData) = c("fibreID", "sampleID", mitochan, channels)
  data = tidyr::pivot_longer(rawData, cols=c(mitochan, channels), names_to="Channel", values_to="Value")

  data_df = as.data.frame(data)

  data_df$sbjType = "control"
  data_df$sbjType[grep("C", data_df$sampleID, invert=TRUE)] = "patient"

  return( data_df )
}











