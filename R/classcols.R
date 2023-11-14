#' @title Classification Colour
#'
#' @description
#' The colour given to a fibre given its probability of being deficient. The
#' colours are on a scale from blue to red where pure blue represents a fibre
#' with a zero probability of being deficient and pure red represents a
#' deficiency probabilty of 1.0.
#'
#' @details
#' The colour gradient is created using [analysis2Dmito::ramp_BlueRed].
#'
#' @param classif A numeric vector of deficiency probabilities, all elements must be between 0.0 and 1.0.
#' @param alphaLevel A number between 0.0 and 1.0 to control the level of the opacity of the colours. Default 1.0.
#'
#' @returns A character vector with hex value RGB names.
#'
#' @importFrom grDevices rgb
#'
#' @export
classcols = function(classif, alphaLevel=1.0){
  alphaCramp = ramp_BlueRed(alphaLevel)
  rgbvals = alphaCramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3], alpha=rgbvals[,4]))
}
