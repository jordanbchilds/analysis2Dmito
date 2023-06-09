#' @title Transparent Colour Creator
#'
#' @description Defining a transparent version of a given colour.
#'
#' @param name A valid colour name, as listed by [colors()].
#' @param rgb A numeic vector of RGB values.
#' @param maxValue The maximum colour value, default to 1.
#' @param alpha The level of transparency
#'
#' @return A character vector of the hexidecimal value defining the colour.
#'
#' @export
#'
#' @importFrom grDevices rgb
#' @importFrom grDevices col2rgb
#' @examples
#' myBlue = alphaCol(0.5, name="blue")
#' myRed = alphaCol(0.1, rgb=c(240, 20, 2-), maxValue=255)
alphaCol = function(alpha, name = NULL, rgbVals = NULL, maxValue=1) {
  if (is.null(name) && is.null(rgbVals)) {
    stop("Either a valid colour name of a vector RGB values must be passed.")
  }
  if (!is.null(name)) {
    rgbVals = as.vector(col2rgb(name))
  }
  return(rgb(rgbVals[1]/maxValue, rgbVals[2]/maxValue, rgbVals[3]/maxValue, alpha))
}
