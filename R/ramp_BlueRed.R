#' @title Blue-Red Colour Ramp
#'
#' @description For creating transitions colours between red and blue, with a given transparency.
#'
#' @details The colour gradient lies between [analysis2Dmito::alphaBlue()] and [analysis2Dmito::alphaRed()]
#'
#' @param alphaLevel The level of opacity of the two end-point colours. A number between 0.0 and 1.0, where 0.0 is completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
ramp_BlueRed = function(alphaLevel){
  colorRamp(c(alphaBlue(alphaLevel), alphaRed(alphaLevel)), alpha=TRUE)
}
