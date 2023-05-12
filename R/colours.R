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
#' @examples
#' myBlue = alphaCol(0.5, name="blue")
#' myRed = alphaCol(0.1, rgb=c(240, 20, 2-), maxValue=255)
alphaCol = function(alpha, name = NULL, rgb = NULL, maxValue=1) {
  if (is.null(name) && is.null(rgb)) {
    stop("Either a valid colour name of a vector RGB values must be passed.")
  }
  if (!is.null(name)) {
    rgb = as.vector(col2rgb(name))
  }
  return(rgb(rgb[1]/maxValue, rgb[2]/maxValue, rgb[3]/maxValue, alpha))
}


#' @title Transparent Black
#'
#' @description Create a black with a specified level of opacity.
#'
#' @details The RGB values for alphaBlack are 0,0,0.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely
#' transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaBlack = function(alpha) rgb(0,0,0, alpha)


#' @title Transparent Blue
#'
#' @description Create a blue with a specified level of opacity.
#'
#' @details The RGB values used for alphaBlue are 0,0 and 128, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaBlue = function(alpha) rgb(0,0,128/255, alpha)


#' @title Transparent Red
#'
#' @description Create a red with a specified level of opacity.
#'
#' @details The RGB values for alphaRed are 255,0 and 0, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaRed = function(alpha) rgb(1,0,0, alpha)


#' @title Transparent Green
#'
#' @description Create a green with a specified level of opacity.
#'
#' @details The RGB values used for alphaGreen are 0,150 and 0, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaGreen = function(alpha) rgb(0,150/255,0, alpha)


#' @title Transparent Pink
#'
#' @description Create a pink with a specified level of opacity.
#'
#' @details The RGB values for alphaPink are 255,62 and 150, where the maximum is 255.
#'
#' @param alpha The level of opacity, between 0.0 and 1.0, 0.0 being completely transparent and 1.0 being opaque.
#'
#' @returns A character object of the hex value colour name created.
#'
#' @export
alphaPink = function(alpha) rgb(255/255,62/255,150/255, alpha)

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


#' @title Classification Colour
#'
#' @description
#' The colour given to a fibre given its probability of being deficient. The
#' colours are on a scale of red to blue where red corresponds to deficiency i.e.
#' a probability of 1.0.
#'
#' @details
#' The colour gradient is created using [analysis2Dmito::ramp_BlueRed()].
#'
#' @param classif A numeric vector for the deficiency probabilities. All elements must be between 0.0 and 1.0.
#' @param alphaLevel A number between 0.0 and 1.0 to control the level of the opacity of the colours.
#'
#' @returns A character vector hex value RGB names.
#'
#' @export
classcols = function(classif, alphaLevel=1.0){
  alphaCramp = ramp_BlueRed(alphaLevel)
  rgbvals = alphaCramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3], alpha=rgbvals[,4]))
}

